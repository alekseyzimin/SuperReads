#ifndef __CHARBUF_HPP__
#define __CHARBUF_HPP__

#include <src/charb.hpp>

/** Streambuf backed by a charb. Similar to std::stringbuf (used by
 * std::ostringstream) but is based on a charb rather than a
 * string. The associated stream (charstream) has an overloaded
 * operator<< which allows to write content since the last rewind.
 */

class charbuf : public std::streambuf {
  charb buf;
public:
  charbuf(size_t s) : std::streambuf(), buf(s) { 
    setp(buf.base(), buf.base() + buf.capacity());
  }
  virtual ~charbuf() { }

  virtual std::streamsize xsputn(const char *s, std::streamsize n) {
    std::streamsize left = epptr() - pptr();
    if(left < n) {
      int off = pptr() - pbase();
      buf.enlarge();
      setp(buf.base(), buf.base() + buf.capacity());
      pbump(off);
    }
    memcpy(pptr(), s, n);
    pbump(n);
    return n;
  }

  virtual int overflow(int c) {
    size_t csize = pptr() - pbase();
    buf.enlarge();
    setp(buf.base(), buf.base() + buf.capacity());
    pbump(csize);
    if(c != EOF) {
      *pptr() = c;
      pbump(1);
    }
    return !EOF;
  }

  size_t size() const { return pptr() - pbase(); }
  char *str() const { return pbase(); }
  void rewind() {
    setp(buf.base(), buf.base() + buf.capacity()); // Any other way to reset pptr?
  }
};

class charstream : public std::ostream {
  charbuf *buf;
public:
  charstream(size_t s = 1024) : std::ostream(new charbuf(s)), buf(0) { 
    buf = (charbuf*)std::ostream::rdbuf();
  }
  virtual ~charstream() {
    std::ostream::rdbuf(0);
    delete buf;
    buf = 0;
  }
  size_t size() const { return buf->size(); }
  const char *str() const { return buf->str(); }
  void rewind() { buf->rewind(); }
};

std::ostream &operator<<(std::ostream &os, const charstream &cs) {
  return os.write(cs.str(), cs.size());
}

#endif /* __CHARBUF_HPP__ */
