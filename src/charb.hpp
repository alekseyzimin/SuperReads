#ifndef __CHARB_HPP__
#define __CHARB_HPP__

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <stdexcept>
#include <string>
#include <cstring>

#include <src/exp_buffer.hpp>

// For the friend business to work
template<typename R> class basic_charb;
template<typename R> char *fgets(basic_charb<R> &b, FILE *stream, char *cptr);
template<typename R> int vsprintf(basic_charb<R> &b, const char *format, va_list ap);
template<typename R> ssize_t getline(basic_charb<R> &b, FILE *stream);
template<typename R> ssize_t getdelim(basic_charb<R> &b, int delim, FILE *stream);

template<typename R>
class basic_charb : public ExpBuffer<char, R> {
  typedef ExpBuffer<char, R> super;
public:
  basic_charb() : super() { }
  explicit basic_charb(size_t s) : super(s + 1) {
    *super::_ptr = '\0';
  }
  basic_charb(const basic_charb &rhs) : super(rhs._base, rhs.capacity()) {
    *super::_ptr = '\0';
  }
  basic_charb(const char *str) : super(str, strlen(str) + 1) {
    *--super::_ptr = '\0';
  }
  basic_charb(const char *str, size_t len) : super(str, len + 1) {
    *--super::_ptr = '\0';
  }
  basic_charb(const std::string &s) : super(s.c_str(), s.size() + 1) {
    *--super::_ptr = '\0';
  }
  virtual ~basic_charb() { }

  basic_charb& operator=(const basic_charb &rhs) {
    super::operator=(rhs);
    *super::_ptr = '\0';
  }

  friend char *fgets <> (basic_charb<R> &b, FILE *stream, char *cptr);
  friend int vsprintf <> (basic_charb<R> &b, const char *format, va_list ap);
  friend ssize_t getline <> (basic_charb<R> &b, FILE *stream);
  friend ssize_t getdelim <> (basic_charb<R> &b, int delim, FILE *stream);
};
typedef basic_charb<reallocator> charb;

/** fgets for char buffer. Expand the size of the buffer if the line
 * does not fit.
 */
template<typename R>
char *fgets(basic_charb<R> &b, FILE *stream, char *cptr) {
  long  pos   = ftell(stream);
  long  npos  = pos;
  char *start = cptr;

  while(true) {
    char *res = fgets(cptr, b.capacity() - (cptr - b._base), stream);
    if(!res)
      break;
    size_t char_read;
    if(pos == -1) {
      char_read = strlen(res);
    } else {
      npos = ftell(stream);
      char_read = npos - pos;
      pos = npos;
    }
    cptr      += char_read;
    if(cptr < b._end - 1 || *(cptr - 1) == '\n')
      break;
    size_t off  = cptr  - b._base;
    size_t soff = start - b._base;
    b.enlarge();
    cptr  = b._base + off;
    start = b._base + soff;
  }
  
  if(cptr == b._base)
    return 0;
  b._ptr = cptr;
  return start;
}

template<typename R>
char *fgets(basic_charb<R> &b, FILE *stream) { return fgets(b, stream, b.base()); }

template<typename R>
char *fgets_append(basic_charb<R> &b, FILE *stream) { return fgets(b, stream, b.ptr()); }

/** Backward compatible fgets for char buffer. The size argument is ignored and present
 * only for backward compatibility.
 */
template<typename T, typename R>
char *fgets(basic_charb<R> &b, T size, FILE *stream) { return fgets(b, stream); }

/** Getline for char buffer.
 */
template<typename R>
ssize_t getline(basic_charb<R> &b, FILE *stream) {
  size_t n = b.capacity();
  ssize_t res = getline(&b._base, &n, stream);
  if(res == -1)
    return res;
  b._ptr = b._base + res;
  if(n != b.capacity())
    b._end = b._base + n;

  return res;
}

template<typename T>
ssize_t getline(basic_charb<T> &b, T *n, int delim, FILE *stream) {
  return getline(b, stream);
}

/** Getdelim for char buffer.
 */
template<typename R>
ssize_t getdelim(basic_charb<R> &b, int delim, FILE *stream) {
  size_t n = b.capacity();
  ssize_t res = getdelim(&b._base, &n, delim, stream);
  if(res == -1)
    return res;
  b._ptr = b._base + res;
  if(n != b.capacity())
    b._end = b._base + n;

  return res;
}

template<typename T, typename R>
ssize_t getdelim(basic_charb<R> &b, T *n, int delim, FILE *stream) {
  return getdelim(b, delim, stream);
}


/** Sprintf for buffer. The buffer grows as needed. The return value
 * is < 0 in case of error, or the number of characters written to the
 * char buffer.
 */
template<typename R>
int sprintf(basic_charb<R> &b, const char *format, ...) {
  va_list ap;

  va_start(ap, format);
  int res = vsprintf(b, format, ap);
  va_end(ap);

  return res;
}

/** Snprintf for backward compatibility.
 */
template<typename T, typename R>
int snprintf(basic_charb<R> &b, T size, const char *format, ...) {
  va_list ap;
  va_start(ap, format);
  int res = vsprintf(b, format, ap);
  va_end(ap);

  return res;
}

/** Vsnprintf for backward compatibility.
 */
template<typename T, typename R>
int vsnprintf(basic_charb<R> &b, T size, const char *format, va_list ap) {
  return vsprintf(b, format, ap);
}

template<typename R>
int vsprintf(basic_charb<R> &b, const char *format, va_list _ap) {
  int res;
  while(true) {
    va_list ap;
    va_copy(ap, _ap);
    res = vsnprintf(b._base, b.capacity(), format, ap);
    va_end(ap);
    if(res < 0)
      return res;
    if((size_t)res < b.capacity())
      break;
    b.ensure(res + 1);
  }
  b._ptr = b._base + res;
  return res;
}
#endif /* __CHARB_HPP__ */
