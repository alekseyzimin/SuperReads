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


#ifndef __REALLOCATORS_HPP__
#define __REALLOCATORS_HPP__

#include <stdlib.h>
#include <sys/mman.h>
#include <string.h>
#include <config.h>

/** Allocators. The interface does not match that of
 * std::allocators. It is composed of only 1 operator which should
 * behave similarly to realloc.
 */

/** Malloc based allocator with no initialization of memory.
 */
template<typename T>
struct reallocator {
  typedef T element_type;
  // T *operator()(T *ptr, size_t osize, size_t nsize) {
  static T* realloc(T *ptr, size_t osize, size_t nsize) {
    return (T*)::realloc((void*)ptr, nsize * sizeof(T));
  }
};

/** Malloc based allocator which also initialized memory.
 */
template<typename T, T v>
struct reallocator_init {
  typedef T element_type;
  static const T default_value = v;
  //  T *operator()(T *ptr, size_t osize, size_t nsize) {
  static T* realloc(T *ptr, size_t osize, size_t nsize) {
    T *res = (T*)::realloc((void*)ptr, nsize * sizeof(T));
    if(res && nsize > osize)
      std::fill(res + osize, res + nsize, v);
    return res;
  }
};
  
#ifdef HAVE_MREMAP
/** Mmap and mremap based allocator. No copying is involved during a
 * realloc. Memory is zeroed out.
 */
template<typename T>
struct remaper {
  typedef T element_type;
  //  T *operator()(T *ptr, size_t osize, size_t nsize) {
  static T* realloc(T *ptr, size_t osize, size_t nsize) {
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
  typedef T element_type;
  static const T default_value = v;
  //  remaper<T> remap;
  //  T *operator()(T *ptr, size_t osize, size_t nsize) {
  static T* realloc(T *ptr, size_t osize, size_t nsize) {
    //    T *res = remap(ptr, osize, nsize);
    T *res = remaper<T>::realloc(ptr, osize, nsize);
    if(res && nsize > osize)
      std::fill((T*)res + osize, (T*)res + nsize, v);
    return res;
  }
};
#else
template<typename T>
struct remaper : public reallocator<T> { 
  typedef T element_type;
  //  T* operator()(T* ptr, size_t osize, size_t nsize) {
  static T* realloc(T* ptr, size_t osize, size_t nsize) {
    T* res = reallocator<T>::operator()(ptr, osize, nsize);
    if(res && nsize > osize)
      memset(res + osize, '\0', sizeof(T) * (nsize - osize));
    return res;
  }
};
template<typename T, T v>
struct remaper_init : public reallocator_init<T, v> {
  typedef T element_type;
};
#endif // HAVE_MREMAP

#endif // __REALLOCATORS_HPP__
