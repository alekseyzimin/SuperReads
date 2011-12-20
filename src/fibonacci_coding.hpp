#ifndef __FIBONACCI_CODING_HPP__
#define __FIBONACCI_CODING_HPP__

#include <errno.h>
#include <stdint.h>
#include <iostream>
#include <algorithm>

struct fibonacci {
  static const uint64_t *fibs;

  template<typename T>
  static int encode(T x, T &res) {
    const uint64_t * f = fibs + 8 * sizeof(T) - 1;
    x++;
    if(x <= (T)0 || x >= *f) {
      errno = EINVAL; // Number is negative or too large
      return -1;
    }
    
    f             = std::upper_bound(fibs, f, x) - 1;
    int high_bit  = f - fibs;
    res           = (T)1 << high_bit;
    x            -= *f;
    while(x > 0) {
      f    = std::upper_bound(fibs, f, x) - 1;
      res |= (T)1 << (f - fibs);
      x   -= *f;
    }
    
    return high_bit + 1;
  }

  template<typename T>
  static T decode(T x) {
    T res = 0;
    for(const uint64_t *f = fibs; x; x >>= 1, ++f)
      if(x & 0x1)
        res += *f;
    return res - 1;
  } 
};

#endif
