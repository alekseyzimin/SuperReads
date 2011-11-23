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
template<typename R> char *strcat(basic_charb<R> &b, const char *src);

template<typename R>
class basic_charb : public ExpBuffer<char, R> {
  typedef ExpBuffer<char, R> super;
public:
  basic_charb() : super() { }
  explicit basic_charb(size_t s) : super(s + 1) {
    *super::ptr_ = '\0';
  }
  basic_charb(const basic_charb &rhs) : super(rhs.base_, rhs.len() + 1) {
    *--super::ptr_ = '\0';
  }
  basic_charb(const char *str) : super(str, strlen(str) + 1) {
    *--super::ptr_ = '\0';
  }
  basic_charb(const char *str, size_t len) : super(str, len + 1) {
    *--super::ptr_ = '\0';
  }
  basic_charb(const std::string &s) : super(s.c_str(), s.size() + 1) {
    *--super::ptr_ = '\0';
  }
  virtual ~basic_charb() { }

  basic_charb &operator=(const char *rhs) {
    size_t rhs_len = strlen(rhs);
    super::reserve(rhs_len + 1);
    strcpy(super::base_, rhs);
    super::ptr_ = super::base_ + rhs_len;
    return *this;
  }

  basic_charb& operator=(basic_charb rhs) {
    this->swap(rhs);
    return *this;
  }
  size_t len() const { return super::size(); }
  friend char *fgets <> (basic_charb<R> &b, FILE *stream, char *cptr);
  friend int vsprintf <> (basic_charb<R> &b, const char *format, va_list ap);
  friend ssize_t getline <> (basic_charb<R> &b, FILE *stream);
  friend ssize_t getdelim <> (basic_charb<R> &b, int delim, FILE *stream);
  friend char *strcat <> (basic_charb<R> &b, const char *src);
};
typedef basic_charb<reallocator<char> > charb;

/** fgets for char buffer. Expand the size of the buffer if the line
 * does not fit.
 */
template<typename R>
char *fgets(basic_charb<R> &b, FILE *stream, char *cptr) {
  long  pos   = ftell(stream);
  long  npos  = pos;
  char *start = cptr;

  while(true) {
    char *res = fgets(cptr, b.capacity() - (cptr - b.base_), stream);
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
    if(cptr < b.end_ - 1 || *(cptr - 1) == '\n')
      break;
    size_t off  = cptr  - b.base_;
    size_t soff = start - b.base_;
    b.enlarge();
    cptr  = b.base_ + off;
    start = b.base_ + soff;
  }
  
  if(cptr == b.base_)
    return 0;
  b.ptr_ = cptr;
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
  ssize_t res = getline(&b.base_, &n, stream);
  if(res == -1)
    return res;
  b.ptr_ = b.base_ + res;
  if(n != b.capacity())
    b.end_ = b.base_ + n;

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
  ssize_t res = getdelim(&b.base_, &n, delim, stream);
  if(res == -1)
    return res;
  b.ptr_ = b.base_ + res;
  if(n != b.capacity())
    b.end_ = b.base_ + n;

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
    res = vsnprintf(b.base_, b.capacity(), format, ap);
    va_end(ap);
    if(res < 0)
      return res;
    if((size_t)res < b.capacity())
      break;
    b.reserve(res + 1);
  }
  b.ptr_ = b.base_ + res;
  return res;
}

/** strcat
 */
template<typename R>
char *strcat(basic_charb<R> &b, const char *src) {
  size_t src_len = strlen(src);
  b.reserve(src_len + b.len() + 1);
  strncpy(b.ptr_, src, src_len + 1);
  b.ptr_ += src_len;
  return b.base_;
}
template<typename T, typename R>
char *strncat(basic_charb<R> &b, T size, const char *src) {
  return strcat(b, src);
}

/** strcpy
 */
template<typename R>
char *strcpy(basic_charb<R> &b, const char *src) {
  b = src;
  return b.begin();
}
template<typename T, typename R>
char *strncpy(basic_charb<R> &b, const char *src) {
  b = src;
  return b.begin();
}
#endif /* __CHARB_HPP__ */
