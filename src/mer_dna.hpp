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
#include <stdint.h>
#include <string.h>

class mer_dna {
public:
  static const uint64_t codes[256];
  static const char     rev_codes[4];
  static const uint64_t bad_code = -1;
    

  // Uninitialized k-mer.
  mer_dna(unsigned int k) : 
    _k(k),
    _data(new uint64_t[k_to_words(k)]) 
  { }

  mer_dna(const mer_dna &m) : 
    _k(m._k),
    _data(new uint64_t[k_to_words(m._k)])
  {
    memcpy(_data, m._data, nb_words() * sizeof(uint64_t));
  }

  // Initialize from a string. The result is unspecified if invalid
  // characters are encountered.
  mer_dna(unsigned int k, const char *s) :
    _k(k),
    _data(new uint64_t[k_to_words(k)])
  {
    from_str(s);
  }

  mer_dna(const std::string &s) :
    _k(s.size()),
    _data(new uint64_t[k_to_words(s.size())])
  {
    from_str(s.c_str());
  }

  ~mer_dna() {
    delete [] _data;
  }

  unsigned int k() const { return _k; }
    
  // Direct access to data. No bound or consistency check. Use with caution!
  //  uint64_t operator[](unsigned int i) { return _data[i]; }
  uint64_t operator[](unsigned int i) const { return _data[i]; }

  bool operator==(const mer_dna& rhs) const;
  bool operator!=(const mer_dna& rhs) { return !this->operator==(rhs); }
  bool operator<(const mer_dna& rhs) const;
  
  class base_proxy {
    uint64_t* const word_;
    unsigned int    i_;
    base_proxy(uint64_t* data, unsigned int i) :
      word_(data + i / (sizeof(uint64_t) * 4)), i_(2 * (i % (sizeof(uint64_t) * 4))) { }
    friend class mer_dna;
  public:
    base_proxy& operator=(char base) { return this->operator=(mer_dna::code(base)); }
    base_proxy& operator=(uint64_t code) {
      uint64_t mask = (uint64_t)0x3 << i_;
      *word_ = (*word_ & ~mask) | (code << i_);
      return *this;
    }
    uint64_t code() const { return (*word_ >> i_) & (uint64_t)0x3; }
    operator char() const { return mer_dna::rev_codes[code()]; }
  };

public:
  base_proxy base(unsigned int i) { return base_proxy(_data, i); }

  // Make current k-mer all As.
  void polyA() {
    memset(_data, '\0', sizeof(uint64_t) * nb_words());
  }

  mer_dna &operator=(const mer_dna &o) {
    if(this != &o) {
      if(_k != o._k)
        throw std::length_error("k-mer of different length");
      memcpy(_data, o._data, sizeof(uint64_t) * nb_words());
    }
    return *this;
  }

  // Shift the k-mer by 1 base, left or right. The char version take
  // a base 'A', 'C', 'G', or 'T'. The uint64_t version takes a code
  // in [0, 3] (not check of validity of argument, taken modulo
  // 4). The return value is the base that was pushed off the side
  // ('N' if the input character is not a valid base).
  uint64_t shift_right(uint64_t c);
  uint64_t shift_left(uint64_t c);

  static uint64_t code(char c) { return codes[(int)c]; }
  static uint64_t complement(uint64_t x) { return (uint64_t)3 - x; }
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
    uint64_t x = codes[(int)c];
    if(x == (uint64_t)-1)
      return 'N';
    return rev_codes[shift_left(x)];
  }

  char shift_right(char c) {
    uint64_t x = codes[(int)c];
    if(x == (uint64_t)-1)
      return 'N';
    return rev_codes[shift_right(x)];
  }

  void reverse_complement();

  // Transform the k-mer into a string. For the char * version,
  // assume that the buffer is large enough to receive k+1
  // characters (space for '\0' at end of string).
  std::string to_str() const;
  void to_str(char *s) const;

  // Get bits [start, start+len). start must be < 2k, len < 64 and
  // start+len < 2k. No checks performed. start and len are in bits, not
  // bases.
  uint64_t get_bits(unsigned int start, unsigned int len) const;


  // Internal stuff

  // Number of words in _data
  unsigned int nb_words() const { return k_to_words(_k); }
  static unsigned int k_to_words(unsigned int k) {
    return (k >> 5) + ((k & 0x1f) != 0);
  }
  // Mask of highest word
  uint64_t msw() const { return ((uint64_t)1 << nb_msb()) - 1; }
  // Nb of bits used in highest word
  unsigned int nb_msb() const {
    uint64_t nb = (_k & 0x1f) << 1;
    return nb ? nb : sizeof(uint64_t) * 8;
  }
  // How much to shift last base in last word of _data
  unsigned int lshift() const { return nb_msb() - 2; }

private:
  static const uint64_t c3 = (uint64_t)0x3;
  const unsigned int    _k;
  uint64_t * const      _data;
    
  bool from_str(const char *s);
  
  // Shift to the right by rs bits (Note bits, not bases)
  void large_shift_right(unsigned int rs);
  static inline uint64_t word_reverse_complement(uint64_t w) {
    w = ((w >> 2)  & 0x3333333333333333UL) | ((w & 0x3333333333333333UL) << 2);
    w = ((w >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((w & 0x0F0F0F0F0F0F0F0FUL) << 4);
    w = ((w >> 8)  & 0x00FF00FF00FF00FFUL) | ((w & 0x00FF00FF00FF00FFUL) << 8);
    w = ((w >> 16) & 0x0000FFFF0000FFFFUL) | ((w & 0x0000FFFF0000FFFFUL) << 16);
    w = ( w >> 32                        ) | ( w                         << 32);
    return ((uint64_t)-1) - w;
  }
};
std::ostream &operator<<(std::ostream &os, const mer_dna &mer);

#endif
