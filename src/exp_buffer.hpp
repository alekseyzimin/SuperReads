#ifndef _EXP_BUFFER_H_
#define _EXP_BUFFER_H_

#include <cstdlib>
#include <stdexcept>
#include <sys/mman.h>
#include <config.h>

/** A growable array (and expanding buffer) but only suitable for
 * types than can be memcpy (e.g., not real constructor). Similar to
 * std::vector but intended to be a drop in replacement of some
 * statically allocated buffers.
 *
 * Potentially the one interest is with very large containers, where
 * it allows growing the size without any construction/destruction of
 * objects (which may not be achievable with
 * std::vector/std::allocator. Is this true?). When using the remaper
 * allocator, no memory copying is even involved.
 *
 * The interface is similar to std::vector by design.
 */


/** Allocators. The interface does not match that of
 * std::allocators. It is composed of only 1 operator which should
 * behave similarly to realloc.
 */

/** Malloc based allocator with no initialization of memory.
 */
template<typename T>
struct reallocator {
  T *operator()(T *ptr, size_t osize, size_t nsize) {
    return (T*)::realloc((void*)ptr, nsize * sizeof(T));
  }
};

/** Malloc based allocator which also initialized memory.
 */
template<typename T, T v>
struct reallocator_init {
  static const T default_value = v;
  T *operator()(T *ptr, size_t osize, size_t nsize) {
    T *res = (T*)::realloc((void*)ptr, nsize * sizeof(T));
    if(res && nsize > osize)
      std::fill(res + osize, res + nsize, v);
  }
};
  
#ifdef HAVE_MREMAP
/** Mmap and mremap based allocator. No copying is involved during a
 * realloc. Memory is zeroed out.
 */
template<typename T>
struct remaper {
  T *operator()(T *ptr, size_t osize, size_t nsize) {
    void *res;
    if(nsize == 0) {     // free
      if(ptr != 0 && osize != 0)
        munmap(ptr, osize * sizeof(T));
      osize = 0;
      return 0;
    }
    if(ptr == 0) {      // alloc
      res = mmap(0, nsize * sizeof(T), PROT_READ|PROT_WRITE, 
                 MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
    } else {            // realloc
      res = mremap(ptr, osize * sizeof(T), nsize * sizeof(T), MREMAP_MAYMOVE);
    }
    if(res == MAP_FAILED)
      return 0;
    return (T*)res;
  }
};

/** Mmap and mremap based allocator with default memory initialization.
 */
template<typename T, T v>
struct remaper_init {
  static const T default_value = v;
  remaper<T> remap;
  T *operator()(T *ptr, size_t osize, size_t nsize) {
    T *res = remap(ptr, osize, nsize);
    if(res && nsize > osize)
      std::fill((T*)res + osize, (T*)res + nsize, v);
    return res;
  }
};
#endif

template<typename T, typename R=reallocator<T> >
class ExpBuffer {
protected:
  T * base_;
  T * end_;
  T * ptr_;
  R realloc;

public:
  typedef T&        reference;
  typedef const T&  const_reference;
  typedef T*        iterator;
  typedef const T*  const_iterator;
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef T         value_type;
  typedef T*        pointer;
  typedef const T*  const_pointer;
  // no reverse iterator type
  // no allocator type (it does not satisfy the std::Allocator interface).

  ExpBuffer() : base_(0), end_(0), ptr_(0) { }
  explicit ExpBuffer(size_type s) : base_(0), end_(0), ptr_(0) {
    reserve(s);
  }
  ExpBuffer(const T *ptr, size_type nb_elements) : base_(0), end_(0), ptr_(0) {
    reserve(nb_elements);
    memcpy(base_, ptr, sizeof(T) * nb_elements);
    ptr_ = base_ + nb_elements;
  }
  ExpBuffer(const ExpBuffer &rhs) : base_(0), end_(0), ptr_(0) {
    reserve(rhs.capacity());
    memcpy(base_, rhs.base_, sizeof(T) * rhs.size());
    ptr_ = base_ + rhs.size();
  }
  virtual ~ExpBuffer() {
    realloc(base_, 0, 0);
  }

  size_type capacity() const { return end_ - base_; }
  size_type size() const { return ptr_ - base_; }

  void swap(ExpBuffer &rhs) {
    std::swap(base_, rhs.base_);
    std::swap(ptr_, rhs.ptr_);
    std::swap(end_, rhs.end_);
  }
  ExpBuffer &operator=(ExpBuffer rhs) {
    swap(rhs);
    return *this;
  }

  template<typename U>
  reference operator[](U i) const { return base_[i]; }
  reference operator*() { return *base_; }
  value_type operator*() const { return *base_; }
  operator T *() const { return base_; }
  bool operator==(const ExpBuffer &rhs) const { return base_ == rhs.base_; }
  bool operator!=(const ExpBuffer &rhs) const { return base_ != rhs.base_; }
  bool operator!() const { return base_ != 0; }
  pointer base() const { return base_; }
  pointer ptr() const { return ptr_; }
  iterator begin() const { return base_; }
  iterator end() const { return ptr_; }
  reference front() const { return *base_; }
  reference back() const { return *(ptr_ - 1); }
  void push_back(const T& x) {
    reserve(size() + 1);
    *ptr_++ = x; // Should we use the inplace new instead?
  }
  void pop_back() {
    if(ptr_ > base_)
      --ptr_;
  }
  void clear() { ptr_ = base_; }
  bool empty() const { return ptr_ == base_; }
  
  void reserve(size_type s) {
    size_type clen = end_ - base_;
    if(s <= clen)
      return;
    if(s <= 2 * clen)
      s = 2 * clen;
    if(s == 0)
      s = 1024;
    T *nbase = (T *)realloc(base_, clen, s);
    if(!nbase)
      throw std::runtime_error("Error allocating memory");
    ptr_  = nbase + (ptr_ - base_);
    end_  = nbase + s;
    base_ = nbase;
  }

  
  void enlarge() { reserve(capacity() * 2); }
};

template<typename T, typename R=reallocator<T> >
class ExpandingBuffer : public ExpBuffer<T, R> {
  typedef ExpBuffer<T, R> super;
  typedef typename super::size_type size_type;
  typedef typename super::pointer pointer;
  typedef typename super::const_pointer const_pointer;
  typedef typename super::reference reference;
public:
  typedef T value_type;
  ExpandingBuffer() : super() { }
  ExpandingBuffer(size_type s) : super(s) { }
  ExpandingBuffer(const ExpandingBuffer &rhs) : super(rhs) { }
  ExpandingBuffer(const_pointer ptr, size_type nb_elements) : 
    super(ptr, nb_elements) { }
  virtual ~ExpandingBuffer() { }

  // Subscript operator which grows the array as needed (similar to
  // Perl's behavior). The new entries are not initialized beyond what
  // the allocator does.
  template<typename U>
  reference operator[](U _i) {
    size_type i = _i;
    if(super::capacity() <= i)
      super::reserve(i + 1);
    if(super::size() <= i)
      super::ptr_ = super::base_ + i + 1;
    return super::operator[](i);
  }
};

// Overloading of the swap operators
namespace std {
  template<typename T, typename R>
  void swap(ExpBuffer<T, R> &a, ExpBuffer<T, R> &b) { a.swap(b); }
  template<typename T, typename R>
  void swap(ExpandingBuffer<T, R> &a, ExpandingBuffer<T, R> &b) { a.swap(b); }
}


#endif /* _EXP_BUFFER_H_ */
