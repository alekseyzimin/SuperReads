/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef _EXP_BUFFER_H_
#define _EXP_BUFFER_H_

#include <assert.h>
#include <cstdlib>
#include <stdexcept>
#include <reallocators.hpp>

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

/* #define CHECK {                      \
    assert(this->base_ <= this->end_); \
    assert(this->base_ <= this->ptr_); \
    assert(this->ptr_ <= this->end_); \
    } while(0); */
#define CHECK

template<typename T, typename R=reallocator<T>, size_t init_size = 1 >
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
    CHECK;
  }
  ExpBuffer(const T *in_ptr, size_type nb_elements) : base_(0), end_(0), ptr_(0) {
    reserve(nb_elements);
    if(in_ptr)
      memcpy(base_, in_ptr, sizeof(T) * nb_elements);
    ptr_ = base_ + nb_elements;
    CHECK;
  }
  ExpBuffer(const ExpBuffer &rhs) : base_(0), end_(0), ptr_(0) {
    reserve(rhs.capacity());
    memcpy(base_, rhs.base_, sizeof(T) * rhs.size());
    ptr_ = base_ + rhs.size();
    CHECK;
  }
  virtual ~ExpBuffer() {
    if(base_)
      realloc(base_, 0, 0);
  }

  void check() const { CHECK; }
  size_type capacity() const { CHECK; return end_ - base_; }
  size_type size() const { CHECK; return ptr_ - base_; }

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
    CHECK;
    if(!ptr_ || ptr_ >= end_) {
      enlarge();
      CHECK;
    }
    *ptr_++ = x; // Should we use the inplace new instead?
    CHECK;
  }
  void pop_back() {
    if(ptr_ > base_)
      --ptr_;
  }
  void clear() { ptr_ = base_; }
  bool empty() const { return ptr_ == base_; }
  bool is_null() const { return !base_ || base_ == end_; }

  void resize(size_type s) {
    reserve(s);
    ptr_ = base_ + s;
  }

  void resize(size_type s, T c) {
    if(s > size()) {
      reserve(s);
      for( ; ptr_ < base_ + s; ++ptr_)
        *ptr_ = c;
    } else
      ptr_ = base_ + s;
  }
  
  void reserve(size_type s = 0) {
    CHECK;
    size_type clen = end_ - base_;
    if(s == 0)
      s = init_size;
    if(s <= clen)
      return;
    if(s <= 2 * clen)
      s = 2 * clen;
    T* nbase = (T*)realloc(base_, clen, s);
    if(!nbase)
      throw std::runtime_error("Error allocating memory");
    ptr_  = nbase + (ptr_ - base_);
    end_  = nbase + s;
    base_ = nbase;
    CHECK;
  }

  
  void enlarge() { reserve(capacity() * 2); }

  /** Touch (read) every page of memory allocated. This forces the
      system to effectively load the pages into memory and can be
      faster when done linearly with 1 thread rather than many threads
      at random.
      The returned value is to be ignored.
   */
  char touch_all() const {
    int  pagesize = getpagesize();
    char ignore   = 0;
    for(const char* ptr = (const char*)base_; ptr < (const char*)end_; ptr += pagesize)
      ignore ^= *ptr;
    return ignore;
  }
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
  ExpandingBuffer(const_pointer in_ptr, size_type nb_elements) : 
    super(in_ptr, nb_elements) { }
  virtual ~ExpandingBuffer() { }

  // Subscript operator which grows the array as needed (similar to
  // Perl's behavior). The new entries are not initialized beyond what
  // the allocator does.
  template<typename U>
  reference operator[](U _i) {
    size_type i = _i;
    assert(super::base_ <= super::end_);
    assert(super::base_ <= super::ptr_ && super::ptr_ <= super::end_);
    if(super::capacity() <= i)
      super::reserve(i + 1);
    assert(super::base_ <= super::end_);
    assert(super::ptr_ <= super::end_);
    if(super::size() <= i)
      super::ptr_ = super::base_ + i + 1;
    assert(super::base_ <= super::end_);
    assert(super::ptr_ <= super::end_);
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
