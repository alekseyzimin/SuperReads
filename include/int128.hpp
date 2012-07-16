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

#ifndef _INT128_H_
#define _INT128_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#ifndef HAVE_INT128
#error "The type __int128 is not supported"
#endif
#endif

// This might be very slow
namespace __int128_ns {
  template<int base>
  void __print_digits(std::ostream& os, unsigned __int128 x,
                      bool lower = true) {
    char buf[50];
    char* ptr = buf + sizeof(buf);
    do {
      int o  = x % base;
      if(o < 10) {
        *--ptr = '0' + o;
      } else {
        *--ptr = (lower ? 'a' : 'A') + o - 10;
      }
      x     /= base;
    } while (x > 0);
    os.write(ptr, buf + sizeof(buf) - ptr);
  }

  template<typename T>
  void __print_decimal(std::ostream& prefix, std::ostream& os, T x,
                       const std::ios::fmtflags& ff) {
    if((ff & std::ios::showpos) && x >= 0)
      prefix << "+";
    if(x == 0) {
      os << "0";
      return;
    }
    if(x < 0) {
      prefix << "-";
      x = -x;
    }
    __print_digits<10>(os, x);
  }

  void __print_bases(std::ostream& prefix, std::ostream& os,
                     unsigned __int128 x, 
                     const std::ios::fmtflags& ff) {
    if(x == 0) {
      os << "0";
      return;
    }
    if(ff & std::ios::showbase) {
      if(ff & std::ios::hex) {
        if(ff & std::ios::uppercase)
          prefix << "0X";
        else
          prefix << "0x";
      } else if(ff & std::ios::oct) {
        prefix << "0";
      }
    }
    if(ff & std::ios::hex) {
      __print_digits<16>(os, (unsigned __int128)x,
                         !(ff & std::ios::uppercase));
    } else if(ff & std::ios::oct) {
      __print_digits<8>(os, (unsigned __int128)x);
    }
  }

  template<typename T>
  void __print_buf(std::ostream& prefix, std::ostream& os, T x,
                   const std::ios::fmtflags& ff) {
    if(ff & std::ios::dec)
      __print_decimal(prefix, os, x, ff);
    else
      __print_bases(prefix, os, (unsigned __int128)x, ff);
  }

  template<typename T>
  void __print(std::ostream&os, T x) {
    const std::ios_base::fmtflags ff = os.flags();

    if(!(ff & std::ios::adjustfield))
      return __print_buf(os, os, x, ff);

    std::ostringstream prefix;
    std::ostringstream buf;
    __print_buf(prefix, buf, x, ff);
    ssize_t nb_padding = os.width() - (prefix.str().size() + buf.str().size());
    if(nb_padding <= 0) {
      os.write(prefix.str().c_str(), prefix.tellp());
      os.write(buf.str().c_str(), buf.tellp());
      return;
    }

    char padding[nb_padding];
    memset(padding, os.fill(), nb_padding);
    if(ff & std::ios::right)
      os.write(padding, nb_padding);
    os.write(prefix.str().c_str(), prefix.tellp());
    if(ff & std::ios::internal)
      os.write(padding, nb_padding);
    os.write(buf.str().c_str(), buf.tellp());
    if(ff & std::ios::left)
      os.write(padding, nb_padding);
  }
}

inline
std::ostream& operator<<(std::ostream& os, __int128 x) {
  __int128_ns::__print(os, x);
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, unsigned __int128 x) {
  __int128_ns::__print(os, x);
  return os;
}

// namespace {
//   template<>
//   class numeric_limits<__int128> {
//   public:
//     static const bool is_specialized = true;
//     static __int128 min() { return (__int128)-1 ^ ((__int128)-1 >> 1); }
//     static __int128 max() { return (__int128)-1 >> 1; }
//     static int digits   = 127;
//     static int digits10 = 38;
    
//   };
// }

#endif /* _INT128_H_ */
