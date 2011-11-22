#ifndef _EXP_BUFFER_H_
#define _EXP_BUFFER_H_

#include <cstdlib>
#include <stdexcept>
#include <sys/mman.h>
#include <config.h>

struct reallocator {
  void *operator()(void *ptr, size_t size) {
    return ::realloc(ptr, size);
  }
};

#ifdef HAVE_MREMAP
struct remaper {
  void *operator()(void *ptr, size_t size) {
    void *res;
    if(size == 0) {     // free
      if(ptr != 0 && csize != 0)
        munmap(ptr, csize);
      csize = 0;
      return 0;
    }
    if(ptr == 0) {      // alloc
      res = mmap(0, size, PROT_READ|PROT_WRITE, 
                 MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
    } else {            // realloc
      res = mremap(ptr, csize, size, MREMAP_MAYMOVE);
    }
    if(res == MAP_FAILED)
      return 0;
    csize = size;
    return res;
  }
  size_t csize;
  remaper() : csize(0) { }
};
#endif

template<typename T, typename R=reallocator>
class ExpBuffer {
protected:
  T * _base;
  T * _end;
  T * _ptr;
  R realloc;

public:
  typedef T   value_type;
  typedef T * iterator;

  ExpBuffer() : _base(0), _end(0), _ptr(0) { }
  explicit ExpBuffer(size_t s) : _base(0), _end(0), _ptr(0) {
    ensure(s);
  }
  ExpBuffer(const T *ptr, size_t nb_elements) : _base(0), _end(0), _ptr(0) {
    ensure(nb_elements);
    memcpy(_base, ptr, sizeof(T) * nb_elements);
    _ptr = _base + nb_elements;
  }
  ExpBuffer(const ExpBuffer &rhs) : _base(0), _end(0), _ptr(0) {
    ensure(rhs.capacity());
    memcpy(_base, rhs._base, sizeof(T) * rhs.size());
    _ptr = _base + rhs.size();
  }
  virtual ~ExpBuffer() {
    if(_base)
      realloc(_base, 0);
    _base = 0;
  }

  size_t capacity() const { return _end - _base; }
  size_t len() const { return _ptr - _base; }

  ExpBuffer &operator=(const ExpBuffer &rhs) {
    if(this != &rhs) {
      ensure(rhs.capacity());
      memcpy(_base, rhs._base, sizeof(T) * rhs.len());
      _ptr = _base + rhs.len();
    }
    return *this;
  }

  template<typename U>
  T& operator[](U i) const { return _base[i]; }
  T& operator*() { return *_base; }
  T operator*() const { return *_base; }
  operator T *() const { return _base; }
  bool operator==(const ExpBuffer &rhs) const { return _base == rhs._base; }
  bool operator!=(const ExpBuffer &rhs) const { return _base != rhs._base; }
  bool operator!() const { return _base != 0; }
  T *base() const { return _base; }
  T *begin() const { return _base; }
  T *ptr() const { return _ptr; }
  T *end() const { return _end; }
  
  void ensure(size_t s) {
    size_t clen = _end - _base;
    if(s <= clen)
      return;
    if(s <= 2 * clen)
      s = 2 * clen;
    if(s == 0)
      s = 1024;
    T *nbase = (T *)realloc(_base, s * sizeof(T));
    if(!nbase)
      throw std::runtime_error("Error allocating memory");
    _ptr  = nbase + (_ptr - _base);
    _end  = nbase + s;
    _base = nbase;
  }
  void enlarge() { ensure(capacity() * 2); }
};

// template<typename T, typename VT>
// class extending_subscript {
//   template<typename U>
//   VT& operator[](U _i) const {
//     T      *self = static_cast<T *>(this);
//     size_t  i    = static_cast<size_t>(_i);
//     if(self->capacity() <= i)
//       self->ensure(i);
//     if(self->len() < i)
//       T::_ptr = T::_ptr + i;
//     return self->base()[i];
//   }
// };

template<typename T, typename R=reallocator>
class ExpandingBuffer : public ExpBuffer<T, R> {
  typedef ExpBuffer<T, R> super;
public:
  typedef T value_type;
  ExpandingBuffer() : super() { }
  ExpandingBuffer(size_t s) : super(s) { }
  ExpandingBuffer(const ExpandingBuffer &rhs) : super(rhs) { }
  ExpandingBuffer(const T *ptr, size_t nb_elements) : super(ptr, nb_elements) { }
  virtual ~ExpandingBuffer() { }

  template<typename U>
  T& operator[](U _i) {
    size_t i = _i;
    if(super::capacity() <= i)
      super::ensure(i + 1);
    if(super::len() <= i)
      super::_ptr = super::_base + i + 1;
    return super::base()[i];
  }
};


#endif /* _EXP_BUFFER_H_ */