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


#include <src/mer_dna.hpp>
#include <iostream>

const uint64_t mer_dna::codes[256] = {
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};
const char mer_dna::rev_codes[4] = { 'A', 'C', 'G', 'T' };

uint64_t mer_dna::shift_left(uint64_t c) {
  const uint64_t     r       = (_data[nb_words()-1] >> lshift()) & c3;
  const unsigned int barrier = nb_words() & ((uint64_t)-1 << 2);
  uint64_t           c2;  // c2 and c1: carries
  uint64_t           c1      = c & c3;
  unsigned int       i       = 0;
  
  for( ; i < barrier; i += 4) {
    c2 = _data[i] >> 62;     _data[i]   = (_data[i] << 2)   | c1;
    c1 = _data[i+1] >> 62;   _data[i+1] = (_data[i+1] << 2) | c2;
    c2 = _data[i+2] >> 62;   _data[i+2] = (_data[i+2] << 2) | c1;
    c1 = _data[i+3] >> 62;   _data[i+3] = (_data[i+3] << 2) | c2;
  }
  c2 = c1;

  switch(nb_words() - i) {
  case 3: c2 = _data[i] >> 62;   _data[i] = (_data[i] << 2) | c1;   ++i;
  case 2: c1 = _data[i] >> 62;   _data[i] = (_data[i] << 2) | c2;   ++i;
  case 1:                        _data[i] = (_data[i] << 2) | c1;
  }
  _data[nb_words() - 1] = _data[nb_words() - 1] & msw();

  return r;
}

uint64_t mer_dna::shift_right(uint64_t c) { 
  const uint64_t r = _data[0] & c3;
  if(nb_words() > 1) {
    const unsigned int barrier = (nb_words() - 1) & ((uint64_t)-1 << 2);
    unsigned int i = 0;
  
    for( ; i < barrier; i += 4) {
      _data[i]   = (_data[i]   >> 2) | (_data[i+1] << 62);
      _data[i+1] = (_data[i+1] >> 2) | (_data[i+2] << 62);
      _data[i+2] = (_data[i+2] >> 2) | (_data[i+3] << 62);
      _data[i+3] = (_data[i+3] >> 2) | (_data[i+4] << 62);
    }
    switch(nb_words() - 1 - i) {
    case 3: _data[i] = (_data[i] >> 2) | (_data[i+1] << 62);  ++i;
    case 2: _data[i] = (_data[i] >> 2) | (_data[i+1] << 62);  ++i;
    case 1: _data[i] = (_data[i] >> 2) | (_data[i+1] << 62);
    }
  }

  _data[nb_words() - 1] = ((_data[nb_words() - 1] & msw()) >> 2) | ((c & c3) << lshift()); 

  return r;
}

void mer_dna::reverse_complement() {
  uint64_t *low  = _data;
  uint64_t *high = _data + nb_words() - 1;
  for( ; low < high; ++low, --high) {
    uint64_t tmp = word_reverse_complement(*low);
    *low         = word_reverse_complement(*high);
    *high        = tmp;
  }
  if(low == high)
    *low = word_reverse_complement(*low);
  unsigned int rs = sizeof(uint64_t) * 8 - nb_msb();
  if(rs > 0)
    large_shift_right(rs);
}

void mer_dna::large_shift_right(unsigned int rs) {
  if(nb_words() > 1) {
    const unsigned int barrier = (nb_words() - 1) & ((uint64_t)-1 << 2);
    const unsigned int ls = sizeof(uint64_t) * 8 - rs;
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
  _data[nb_words() - 1] = (_data[nb_words() - 1] >> rs) & msw();
}

std::string mer_dna::to_str() const {
  char s[_k + 1];
  to_str(s);
  return std::string(s, _k);
}

void mer_dna::to_str(char *s) const {
  char *cs = s + _k;
  *cs-- = '\0';
  unsigned int j;
  for(j = 0; j < nb_words() - 1; ++j) {
    uint64_t x = _data[j];
    for(unsigned int k = 0; k < 32; ++k, --cs, x >>= 2)
      *cs = rev_codes[x & c3];
  }
  uint64_t x = _data[j];
  for(unsigned int k = 0; k < (_k & 0x1f); ++k, --cs, x >>= 2)
    *cs = rev_codes[x & c3];
}

// static uint64_t reverse_complement_word(uint64_t v, uint_t length) {
//   v = ((v >> 2)  & 0x3333333333333333UL) | ((v & 0x3333333333333333UL) << 2);
//   v = ((v >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((v & 0x0F0F0F0F0F0F0F0FUL) << 4);
//   v = ((v >> 8)  & 0x00FF00FF00FF00FFUL) | ((v & 0x00FF00FF00FF00FFUL) << 8);
//   v = ((v >> 16) & 0x0000FFFF0000FFFFUL) | ((v & 0x0000FFFF0000FFFFUL) << 16);
//   v = ( v >> 32                        ) | ( v                         << 32);
//   return (((uint64_t)-1) - v) >> ((8*sizeof(v)) - (length << 1));
// }

// void mer_dna::reverse_complement() {
//   uint64_t *right = _data;
//   uint64_t *left = _data + (nb_words() - 1);

//   while(left > right) {
    
//   }
// }

bool mer_dna::from_str(const char *s) {
  const char *cs = s + _k - 1;
  for(unsigned int j = 0; j < nb_words(); ++j) {
    uint64_t *x = &_data[j];
    *x = 0;
    for(unsigned int k = 0; cs >= s && k < 32; ++k, --cs) {
        uint64_t c = codes[(int)*cs];
      if(c == (uint64_t)-1)
        return false;
        *x |= (c << (2*k));
    }
  }
  return true;
}

std::ostream &operator<<(std::ostream &os, const mer_dna &mer) {
  char s[mer.k() + 1];
  mer.to_str(s);
  return os << s;
}

uint64_t mer_dna::get_bits(unsigned int start, unsigned int len) const {
  uint64_t res = _data[start >> 6] >> (start & 0x3f);
  if(len > 64 - (start & 0x3f))
    res |= _data[(start >> 6) + 1] << (64 - (start & 0x3f));
  return res & (((uint64_t)1 << len) - 1);
}

bool mer_dna::operator==(const mer_dna& rhs) const  {
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

bool mer_dna::operator<(const mer_dna& rhs) const {
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
