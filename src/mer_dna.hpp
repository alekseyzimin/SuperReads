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


#ifndef __JELLYFISH_MER_DNA_HPP__
#define __JELLYFISH_MER_DNA_HPP__

#include <iostream>
#include <string>
#include <stdexcept>
#include <limits>
#include <stdint.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace mer_dna_ns {
  struct dna_codes {
    static const int  codes[256];
    static const char rev_codes[4];
  };


  // Checkered mask. cmask<uint16_t, 1> is every other bit on
  // (0x55). cmask<uint16_t,2> is two bits one, two bits off (0x33). Etc.
  template<typename U, int len, int l = sizeof(U) * 8 / (2 * len)>
  struct cmask {
    static const U v =
      (cmask<U, len, l - 1>::v << (2 * len)) | (((U)1 << len) - 1);
  };
  template<typename U, int len>
  struct cmask<U, len, 0> {
    static const U v = 0;
  };

  // Fast reverse complement of one word through bit tweedling.
  inline uint64_t word_reverse_complement(uint32_t w) {
    typedef uint64_t U;
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ( w >> 16                   ) | ( w                    << 16);
    return ((U)-1) - w;
  }

  inline uint64_t word_reverse_complement(uint64_t w) {
    typedef uint64_t U;
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
    w = ( w >> 32                   ) | ( w                    << 32);
    return ((U)-1) - w;
  }

#ifdef HAVE_INT128
  inline unsigned __int128 word_reverse_complement(unsigned __int128 w) {
    typedef unsigned __int128 U;
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
    w = ((w >> 32) & cmask<U, 32>::v) | ((w & cmask<U, 32>::v) << 32);
    w = ( w >> 64                   ) | ( w                    << 64);
    return ((U)-1) - w;
  }
#endif

  template <typename T = uint64_t>
  class mer_base {
  public:
    enum { CODE_A, CODE_C, CODE_G, CODE_T,
           CODE_RESET = -1, CODE_IGNORE = -2, CODE_COMMENT = -3 };
    typedef T base_type;
    
    // Uninitialized k-mer.
    mer_base(unsigned int k) : 
      _k(k),
      _data(new base_type[nb_words()]) 
    { }

    mer_base(const mer_base &m) : 
      _k(m.k()),
      _data(new base_type[nb_words()])
    {
      memcpy(_data, m._data, nb_words() * sizeof(base_type));
    }

    // Initialize from a string. The result is unspecified if invalid
    // characters are encountered.
    mer_base(unsigned int k, const char *s) :
      _k(k),
      _data(new base_type[nb_words()])
    {
      from_chars(s);
    }

    mer_base(const std::string &s) :
      _k(s.size()),
      _data(new base_type[nb_words()])
    {
      from_chars(s.begin());
    }

    ~mer_base() {
      delete [] _data;
    }

    unsigned int k() const { return _k; }
    
    // Direct access to data. No bound or consistency check. Use with caution!
    //  base_type operator[](unsigned int i) { return _data[i]; }
    base_type word(unsigned int i) const { return _data[i]; }
    base_type operator[](unsigned int i) const { return _data[i]; }

    bool operator==(const mer_base& rhs) const {
      if(_k != rhs._k)
        return false;
      unsigned int i = nb_words() - 1;
      bool res = (_data[i] & msw()) == (rhs._data[i] & msw());
      while(res && i > 7) {
        i -= 8;
        res &= _data[i+7] == rhs._data[i+7]; 
        res &= _data[i+6] == rhs._data[i+6]; 
        res &= _data[i+5] == rhs._data[i+5]; 
        res &= _data[i+4] == rhs._data[i+4]; 
        res &= _data[i+3] == rhs._data[i+3]; 
        res &= _data[i+2] == rhs._data[i+2];
        res &= _data[i+1] == rhs._data[i+1];
        res &= _data[i]   == rhs._data[i];
      }
      switch(i) {
      case 7: res &= _data[6] == rhs._data[6];
      case 6: res &= _data[5] == rhs._data[5];
      case 5: res &= _data[4] == rhs._data[4];
      case 4: res &= _data[3] == rhs._data[3];
      case 3: res &= _data[2] == rhs._data[2];
      case 2: res &= _data[1] == rhs._data[1];
      case 1: res &= _data[0] == rhs._data[0];
      }
      return res;
    }

    bool operator!=(const mer_base& rhs) { return !this->operator==(rhs); }
    bool operator<(const mer_base& rhs) const {
      if(_k != rhs._k)
        return false;  // Really they are not comparable. Should we throw instead?
      unsigned int i = nb_words();
      while(i >= 8) {
        i -= 8;
        if(_data[i+7] != rhs._data[i+7]) return _data[i+7] < rhs._data[i+7];
        if(_data[i+6] != rhs._data[i+6]) return _data[i+6] < rhs._data[i+6];
        if(_data[i+5] != rhs._data[i+5]) return _data[i+5] < rhs._data[i+5];
        if(_data[i+4] != rhs._data[i+4]) return _data[i+4] < rhs._data[i+4];
        if(_data[i+3] != rhs._data[i+3]) return _data[i+3] < rhs._data[i+3];
        if(_data[i+2] != rhs._data[i+2]) return _data[i+2] < rhs._data[i+2];
        if(_data[i+1] != rhs._data[i+1]) return _data[i+1] < rhs._data[i+1];
        if(_data[i]   != rhs._data[i])   return _data[i]   < rhs._data[i];
      }
      switch(i) {
      case 7: if(_data[6] != rhs._data[6]) return _data[6] < rhs._data[6];
      case 6: if(_data[5] != rhs._data[5]) return _data[5] < rhs._data[5];
      case 5: if(_data[4] != rhs._data[4]) return _data[4] < rhs._data[4];
      case 4: if(_data[3] != rhs._data[3]) return _data[3] < rhs._data[3];
      case 3: if(_data[2] != rhs._data[2]) return _data[2] < rhs._data[2];
      case 2: if(_data[1] != rhs._data[1]) return _data[1] < rhs._data[1];
      case 1: if(_data[0] != rhs._data[0]) return _data[0] < rhs._data[0];
      }
      return false;
    }

  
    class base_proxy {
      base_type* const word_;
      unsigned int    i_;
      base_proxy(base_type* data, unsigned int i) :
        word_(data + i / wbases), i_(2 * (i % wbases)) { }
      friend class mer_base;
    public:
      base_proxy& operator=(char base) { return this->operator=(mer_base::code(base)); }
      base_proxy& operator=(int code) {
        base_type mask = c3 << i_;
        *word_ = (*word_ & ~mask) | ((base_type)code << i_);
        return *this;
      }
      int code() const { return (*word_ >> i_) & (base_type)0x3; }
      operator char() const { return dna_codes::rev_codes[code()]; }
    };

  public:
    base_proxy base(unsigned int i) { return base_proxy(_data, i); }

    // Make current k-mer all As.
    void polyA() {
      memset(_data, '\0', sizeof(base_type) * nb_words());
    }

    mer_base &operator=(const mer_base &o) {
      if(this != &o) {
        if(_k != o._k)
          throw std::length_error("k-mer of different length");
        memcpy(_data, o._data, sizeof(base_type) * nb_words());
      }
      return *this;
    }

    mer_base& operator=(const char* s) {
      if(strlen(s) < _k)
        throw std::length_error("Input string is to short");
      from_chars(s);
      return *this;
    }

    mer_base& operator=(const std::string& s) {
      if(s.size() < _k)
        throw std::length_error("Input string is to short");
      from_chars(s.c_str());
      return *this;
    }

    // Shift the k-mer by 1 base, left or right. The char version take
    // a base 'A', 'C', 'G', or 'T'. The base_type version takes a code
    // in [0, 3] (not check of validity of argument, taken modulo
    // 4). The return value is the base that was pushed off the side
    // ('N' if the input character is not a valid base).
    base_type shift_left(int c) {
      const base_type    r       = (_data[nb_words()-1] >> lshift()) & c3;
      const unsigned int barrier = nb_words() & (~c3);
      base_type          c2;    // c2 and c1: carries
      base_type          c1      = (base_type)c & c3;
      unsigned int       i       = 0;
  
      for( ; i < barrier; i += 4) {
        c2 = _data[i]   >> wshift;   _data[i]   = (_data[i]   << 2) | c1;
        c1 = _data[i+1] >> wshift;   _data[i+1] = (_data[i+1] << 2) | c2;
        c2 = _data[i+2] >> wshift;   _data[i+2] = (_data[i+2] << 2) | c1;
        c1 = _data[i+3] >> wshift;   _data[i+3] = (_data[i+3] << 2) | c2;
      }
      c2 = c1;

      switch(nb_words() - i) {
      case 3: c2 = _data[i] >> wshift;   _data[i] = (_data[i] << 2) | c1;   ++i;
      case 2: c1 = _data[i] >> wshift;   _data[i] = (_data[i] << 2) | c2;   ++i;
      case 1:                        _data[i] = (_data[i] << 2) | c1;
      }
      _data[nb_words() - 1] &= msw();

      return r;
    }

    base_type shift_right(int c) { 
      const base_type r = _data[0] & c3;
      if(nb_words() > 1) {
        const unsigned int barrier = (nb_words() - 1) & (~c3);
        unsigned int i = 0;
  
        for( ; i < barrier; i += 4) {
          _data[i]   = (_data[i]   >> 2) | (_data[i+1] << wshift);
          _data[i+1] = (_data[i+1] >> 2) | (_data[i+2] << wshift);
          _data[i+2] = (_data[i+2] >> 2) | (_data[i+3] << wshift);
          _data[i+3] = (_data[i+3] >> 2) | (_data[i+4] << wshift);
        }
        switch(nb_words() - 1 - i) {
        case 3: _data[i] = (_data[i] >> 2) | (_data[i+1] << wshift);  ++i;
        case 2: _data[i] = (_data[i] >> 2) | (_data[i+1] << wshift);  ++i;
        case 1: _data[i] = (_data[i] >> 2) | (_data[i+1] << wshift);
        }
      }

      _data[nb_words() - 1] = 
        ((_data[nb_words() - 1] & msw()) >> 2) | (((base_type)c & c3) << lshift()); 

      return r;
    }

    // Non DNA codes are negative
    static inline bool not_dna(int c) { return c < 0; }
    static int code(char c) { return dna_codes::codes[(int)(unsigned char)c]; }
    static char rev_code(int x) { return dna_codes::rev_codes[x]; }
    static int complement(int x) { return (base_type)3 - x; }
    static char complement(char c) { 
      switch(c) {
      case 'A': case 'a': return 'T'; 
      case 'C': case 'c': return 'G';
      case 'G': case 'g': return 'C';
      case 'T': case 't': return 'A';
      }
      return 'N';
    }

    char shift_left(char c) { 
      int x = code(c);
      if(x == -1)
        return 'N';
      return rev_code(shift_left(x));
    }

    char shift_right(char c) {
      int x = code(c);
      if(x == -1)
        return 'N';
      return rev_code(shift_right(x));
    }

    void reverse_complement() {
      base_type *low  = _data;
      base_type *high = _data + nb_words() - 1;
      for( ; low < high; ++low, --high) {
        base_type tmp = word_reverse_complement(*low);
        *low          = word_reverse_complement(*high);
        *high         = tmp;
      }
      if(low == high)
        *low = word_reverse_complement(*low);
      unsigned int rs = wbits - nb_msb();
      if(rs > 0)
        large_shift_right(rs);
    }

    void canonicalize() {
      mer_base rc = this->get_reverse_complement();
      if(rc < *this)
        *this = rc;
    }

    mer_base get_reverse_complement() const {
      mer_base res(*this);
      res.reverse_complement();
      return res;
    }

    mer_base get_canonical() const {
      mer_base rc = this->get_reverse_complement();
      if(rc < *this)
        return rc;
      else
        return *this;
    }

    // Transfomr the k-mer into a C++ string.
    std::string to_str() const {
      std::string res(_k, '\0');
      to_chars(res.begin());
      return res;
    }

    // Transform the k-mer into a string. For the char * version,
    // assume that the buffer is large enough to receive k+1
    // characters (space for '\0' at end of string).
    void to_str(char* s) const {
      s = to_chars(s);
      *s = '\0';
    }

    // Copy bases as char to the output iterator it. No '\0' is added
    // or check made that there is enough space. The iterator pointed
    // after the last base is returned.
    template<typename OutputIterator>
    OutputIterator to_chars(OutputIterator it) const {
      int shift  = lshift(); // Number of bits to shift to get base
      
      for(int j = nb_words() - 1; j >= 0; --j) {
        base_type w = _data[j];
        for( ; shift >= 0; shift -= 2, ++it)
          *it = rev_code((w >> shift) & c3);
        shift = wshift;
      }
      return it;
    }

    // Get bits [start, start+len). start must be < 2k, len <
    // sizeof(base_type) and start+len < 2k. No checks
    // performed. start and len are in bits, not bases.
    base_type get_bits(unsigned int start, unsigned int len) const {
      unsigned int q = start / wbits;
      unsigned int r = start % wbits;

      base_type res = _data[q] >> r;
      if(len > wbits - r)
        res |= _data[q + 1] << (wbits - r);
      return res & (((base_type)1 << len) - 1);
    }


    // Internal stuff

    // Number of words in _data
    inline unsigned int nb_words() const { 
      return (_k / wbases) + (_k % wbases != 0);
    }

    // Mask of highest word
    inline base_type msw() const { return ((base_type)1 << nb_msb()) - 1; }
    // Nb of bits used in highest word
    inline  unsigned int nb_msb() const {
      base_type nb = (_k % wbases) * 2;
      return nb ? nb : wbits;
    }
    // How much to shift last base in last word of _data
    inline unsigned int lshift() const { return nb_msb() - 2; }

  private:
    static const base_type c3 = (base_type)0x3;
    static const int wshift = sizeof(base_type) * 8 - 2; // left shift in 1 word
    static const int wbases = 4 * sizeof(base_type); // bases in a word
    static const int wbits  = 8 * sizeof(base_type); // bits in a word
    const unsigned int     _k;
    base_type * const      _data;
    
    template<typename InputIterator>
    bool from_chars(InputIterator it) {
      int shift = lshift();

      for(int j = nb_words() - 1; j >= 0; --j) {
        base_type& w = _data[j];
        w = 0;
        for( ; shift >= 0; shift -= 2, ++it) {
          int c = code(*it);
          if(not_dna(c))
            return false;
          w |= (base_type)c << shift;
        }
        shift = wshift;
      }
      return true;
    }

    // Shift to the right by rs bits (Note bits, not bases)
    void large_shift_right(unsigned int rs) {
      if(nb_words() > 1) {
        const unsigned int barrier = (nb_words() - 1) & (~c3);
        const unsigned int ls = wbits - rs;
        unsigned int i = 0;

        for( ; i < barrier; i += 4) {
          _data[i]   = (_data[i]   >> rs) | (_data[i+1] << ls);
          _data[i+1] = (_data[i+1] >> rs) | (_data[i+2] << ls);
          _data[i+2] = (_data[i+2] >> rs) | (_data[i+3] << ls);
          _data[i+3] = (_data[i+3] >> rs) | (_data[i+4] << ls);
        }
        switch(nb_words() - 1 - i) {
        case 3: _data[i] = (_data[i] >> rs) | (_data[i+1] << ls); ++i;
        case 2: _data[i] = (_data[i] >> rs) | (_data[i+1] << ls); ++i;
        case 1: _data[i] = (_data[i] >> rs) | (_data[i+1] << ls);
        }
      }
      _data[nb_words() - 1] >>= rs;
      _data[nb_words() - 1]  &= msw();
    }
  };
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const mer_dna_ns::mer_base<T>& mer) {
  char s[mer.k() + 1];
  mer.to_str(s);
  return os << s;
}

typedef mer_dna_ns::mer_base<uint32_t> mer_dna32;
typedef mer_dna_ns::mer_base<uint64_t> mer_dna64;
#ifdef HAVE_INT128
typedef mer_dna_ns::mer_base<unsigned __int128> mer_dna128;
#endif

typedef mer_dna64 mer_dna;


#endif
