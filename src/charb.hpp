#ifndef __CHARB_HPP__
#define __CHARB_HPP__

#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <string>

class charb {
  char *_base;
  char *_end;
  char *_ptr;

public:
  charb() : _base(0), _end(0), _ptr(0) { }
  explicit charb(size_t s) : _base(0), _end(0), _ptr(0) {
    ensure(s);
    *_base = '\0';
  }
  charb(const charb &rhs) : _base(0), _end(0), _ptr(0) {
    ensure(rhs.capacity());
    memcpy(_base, rhs._base, capacity());
  }
  charb(const char *str) : _base(0), _end(0), _ptr(0) {
    size_t len = strlen(str);
    ensure(len + 1);
    memcpy(_base, str, len + 1);
    _ptr = _base + len;
  }
  charb(const char *str, size_t len) : _base(0), _end(0), _ptr(0) {
    ensure(len + 1);
    memcpy(_base, str, len + 1);
    _ptr = _base + len;
  }
  charb(const std::string &s) : _base(0), _end(0), _ptr(0) {
    ensure(s.size() + 1);
    memcpy(_base, s.c_str(), s.size() + 1);
    _ptr = _base + s.size();
  }
  virtual ~charb() {
    if(_base)
      free(_base);
  }

  size_t capacity() const { return _end - _base; }
  size_t len() const { return _ptr - _base; }

  charb& operator=(const charb &rhs) {
    if(this == &rhs)
      return *this;
    ensure(rhs.capacity());
    memcpy(_base, rhs._base, capacity());
  }

  char& operator[](size_t i) { return _base[i]; }
  char operator[](size_t i) const { return _base[i]; }
  char& operator*() { return *_base; }
  char operator*() const { return *_base; }
  operator char *() const { return _base; }
  operator const char *() const { return _base; }

  void ensure(size_t s) {
    size_t clen = _end - _base;
    if(s <= clen)
      return;
    if(s <= 2 * clen)
      s = 2 * clen;
    char *nbase = (char *)realloc(_base, s);
    if(!nbase)
      throw std::runtime_error("Error allocating memory");
    _ptr  = nbase + (_ptr - _base);
    _end  = nbase + s;
    _base = nbase;
  }
  void enlarge() { ensure(capacity() * 2); }

  friend char *fgets(charb &b, FILE *stream);
};

/** fgets for char buffer. Expand the size of the buffer if the line
 * does not fit.
 */
char *fgets(charb &b, FILE *stream) {
  char *cptr = b._base;
  long  pos  = ftell(stream);
  
  while(true) {
    char *res = fgets(cptr, b.capacity() - (cptr - b._base), stream);
    if(!res)
      break;
    long npos = ftell(stream);
    cptr      += npos - pos;
    if(cptr < b._end - 1 || *(cptr - 1) == '\n')
      break;
    size_t off = cptr - b._base;
    b.enlarge();
    cptr = b._base + off;
    pos = npos;
  }
  
  if(cptr == b._base)
    return 0;
  b._ptr = cptr;
  return b._base;
}

/** Backward compatible fgets for char buffer. The size argument is ignored and present
 * only for backward compatibility.
 */
char *fgets(charb &b, int size, FILE *stream) { return fgets(b, stream); }

#endif /* __CHARB_HPP__ */
