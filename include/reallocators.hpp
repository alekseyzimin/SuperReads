#ifndef __REALLOCATORS_HPP__
#define __REALLOCATORS_HPP__

#include <sys/mman.h>
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
  T *operator()(T *ptr, size_t osize, size_t nsize) {
    return (T*)::realloc((void*)ptr, nsize * sizeof(T));
  }
};

/** Malloc based allocator which also initialized memory.
 */
template<typename T, T v>
struct reallocator_init {
  typedef T element_type;
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
  typedef T element_type;
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
  typedef T element_type;
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

#endif // __REALLOCATORS_HPP__
