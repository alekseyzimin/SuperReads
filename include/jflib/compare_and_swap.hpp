#ifndef _JFLIB_COMPARE_AND_SWAP_H_
#define _JFLIB_COMPARE_AND_SWAP_H_

#include <sys/types.h>
#include <stdint.h>
#include <jflib/atomic_field.hpp>

namespace jflib {
  // Numeric type of length rounded up to the size of a
  // word. Undefined, and raise a compilation error, if the length is
  // not a machine word size
  template<typename T, size_t n> union word_t;
  template<typename T> union word_t<T, 1> { typedef uint8_t  w_t; T v; w_t w; };
  template<typename T> union word_t<T, 2> { typedef uint16_t w_t; T v; w_t w; };
  template<typename T> union word_t<T, 4> { typedef uint32_t w_t; T v; w_t w; };
  template<typename T> union word_t<T, 8> { typedef uint64_t w_t; T v; w_t w; };
  
  /** Type safe version of CAS.
   * @param [in] ptr Memory location.
   * @param [in] ov  Presumed value at location.
   * @param [in] nv  Value to write.
   * @param [out] cv Value at location at time of call.
   * @return true if CAS is successful.
   *
   * The CAS operation is successful if, at the time of call, ov is
   * equal to *ptr, the value at the memory location. In that case, nv
   * is written to *ptr, and when the call returns, cv == ov.
   *
   * If it fails, cv contains *ptr at the time of call.
   */
  template<typename T>
  bool cas(T *ptr, const T &ov, const T &nv, T *cv) {
    typedef word_t<T, sizeof(T)> val_t;
    val_t _cv, _ov, _nv;
    _ov.v = ov;
    _nv.v = nv;
    _cv.w = __sync_val_compare_and_swap((typename val_t::w_t *)ptr, _ov.w, _nv.w);
    *cv = _cv.v;
    return _cv.w == _ov.w;
  }

  /** Type safe version of CAS.  Identical to 4 argument version,
   * except does not return the previous value.
   */
  template<typename T>
  bool cas(T *ptr, const T &ov, const T &nv) {
    T cv;
    return cas(ptr, ov, nv, &cv);
  }
}


#endif /* _JFLIB_COMPARE_AND_SWAP_H_ */
