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

#ifndef __DIVISOR_HPP__
#define __DIVISOR_HPP__

#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/// Fast division and modulo by a fix positive number. This class speeds up
/// repeated computation of division and modulo by the same
/// number. For example, if a hash table, for every operation one
/// computes `hash(value) % size`. The computation of the modulo by
/// this constant size can be improved by this class. @see divisor64.
///
/// This class requires the type __int128 (define HAVE_INT128). If not
/// supported, it will use the standard (slow) division and modulo
/// operator.
///
/// Example:
/// ~~~~{.cc}
/// divisor64 d(7) // Fast divisor by 7
/// uint64_t r = 0;
/// for(uint64_t i = 0; i < 1000000; ++i)
///   r += i % d;
/// std::cout << r << "\n";
/// ~~~~

namespace divisor {
inline int leading_zeroes(int x) { return __builtin_clz(x); } // CLK
inline int leading_zeroes(unsigned int x) { return __builtin_clz(x); }
inline int leading_zeroes(unsigned long x) { return __builtin_clzl(x); }
inline int leading_zeroes(unsigned long long x) { return __builtin_clzll(x); }

template <typename T>
unsigned int floorLog2(T n) {
  return sizeof(T) * 8 - 1 - leading_zeroes(n);
}

template<typename T>
unsigned int ceilLog2(T n) {
  unsigned int r = floorLog2(n);
  return n > (((T)1) << r) ? r + 1 : r;
}

template<typename T>
T div_ceil(T a, T b) {
  T q = a / b;
  return a % b == 0 ? q : q + 1;
}

/// Fast divisor class.
class divisor64 {
  const uint64_t d_;
#ifdef HAVE_INT128
  const uint16_t p_;
  const unsigned __int128 m_;
#endif

public:
  /// Constructor to divide by `d` or get modulo `d`.
  ///
  /// @param d The divisor
  divisor64(uint64_t d) :
    d_(d)
#ifdef HAVE_INT128
    , p_(ceilLog2(d_)),
    m_((div_ceil((unsigned __int128)1 << (64 + p_), (unsigned __int128)d_)) & (uint64_t)-1)
#endif
  { }

  /// Divide `n` by `d`.
  ///
  /// @param n The number to divide.
  /// @return `n / d`
  inline uint64_t divide(const uint64_t n) const {
#ifdef HAVE_INT128
    switch(m_) {
    case 0:
      return n >> p_;
    default:
      const unsigned __int128 n_ = (unsigned __int128)n;
      return (n_ + ((n_ * m_) >> 64)) >> p_;
    }
#else
    return n / d_;
#endif
  }

  /// Get `n` modulo `d`.
  ///
  /// @param n The number to get the modulo of.
  /// @return `n % d`
  inline uint64_t remainder(uint64_t n) const {
#ifdef HAVE_INT128
    switch(m_) {
    case 0:
      return n & (((uint64_t)1 << p_) - 1);
    default:
      return n - divide(n) * d_;
    }
#else
    return n % d_;
#endif
  }

  /// Euclidian division. Sets `q <- n / d` and `r <- n % d`. This is
  /// faster than doing each independently.
  ///
  /// @param [in] n The number to do the Euclidian division with `d`
  /// @param [out] q Quotient
  /// @param [out] r Rest
  inline void division(uint64_t n, uint64_t &q, uint64_t &r) const {
#ifdef HAVE_INT128
    switch(m_) {
    case 0:
      q = n >> p_;
      r = n & (((uint64_t)1 << p_) - 1);
      break;
    default:
      q = divide(n);
      r = n - q * d_;
      break;
    }
#else
    q = n / d_;
    r = n % d_;
#endif
  }

  /// Divisor `d`
  uint64_t d() const { return d_; }
  /// Internal shift value `p`
  uint64_t p() const {
#ifdef HAVE_INT128
    return p_;
#else
    return 0;
#endif
  }
  /// Internal multipler value `m`
  uint64_t m() const {
#ifdef HAVE_INT128
    return m_;
#else
    return 0;
#endif
  }
};

/// Division operator. Divides `n` by `d`
inline uint64_t operator/(uint64_t n, const divisor64& d) {
  return d.divide(n);
}
/// Modulo operator. Computes `n % d`
inline uint64_t operator%(uint64_t n, const divisor64& d) {
  return d.remainder(n);
}
}

typedef divisor::divisor64 divisor64;

inline std::ostream& operator<<(std::ostream& os, const divisor64& d) {
  return os << "d:" << d.d() << ",p:" << d.p() << ",m:" << d.m();
}

#endif /* __DIVISOR_HPP__ */
