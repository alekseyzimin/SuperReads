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

#ifndef __DNA_FRAGMENT_HPP__
#define __DNA_FRAGMENT_HPP__

#include <jellyfish/mer_dna.hpp>

// DNA fragment. Unlike a mer, the length can change during the life
// time of the object.
template<typename T = uint64_t>
class dna_fragment_base :
  public jellyfish::mer_dna_ns::mer_base<T, dna_fragment_base<T> >
{
public:
  typedef typename jellyfish::mer_dna_ns::mer_base<T, dna_fragment_base<T> > super;
  typedef T base_type;

  dna_fragment_base() : super(0), k_(0) { }
  explicit dna_fragment_base(unsigned int k) : super(k), k_(k) { }
  dna_fragment_base(const dna_fragment_base& rhs) : super(rhs), k_(rhs.k()) { }
  dna_fragment_base(dna_fragment_base&& rhs) : super(0), k_(rhs.k()) {
    std::swap(super::_data, rhs._data);
  }
  dna_fragment_base(unsigned int k, const char* s) : super(k), k_(k) {
    super::from_chars(s);
  }
  explicit dna_fragment_base(const char* s) : super(strlen(s)), k_(strlen(s)) {
    super::from_chars(s);
  }
  explicit dna_fragment_base(const std::string& s) : super(s.size()), k_(s.size()) {
    super::from_chars(s.begin());
  }
  // Raw representation
  dna_fragment_base(unsigned int k, void* raw) : super(k, raw), k_(k) { }

  ~dna_fragment_base() { }

  dna_fragment_base& operator=(const dna_fragment_base& rhs) {
    super::_data = super::allocator::realloc(super::_data, super::nb_words(k_),
                                             super::nb_words(rhs.k_));
    k_ = rhs.k_;
    super::operator=(rhs);
    return *this;
  }
  dna_fragment_base& operator=(dna_fragment_base&& rhs) {
    std::swap(super::_data, rhs._data);
    k_ = rhs.k_;
    return *this;
  }

  /// Get the k-mer starting at a given base of the fragment. If start
  /// + k is greater than the length of the fragment, then the output
  /// k-mer is only partially (or not at all) filled. The part of the
  /// k-mer which extends beyond the fragment is random.
  template<typename U = T>
  jellyfish::mer_dna_ns::mer_base_static<U> sub_mer(const unsigned int start) const {
    jellyfish::mer_dna_ns::mer_base_static<U> res;
    const unsigned k = jellyfish::mer_dna_ns::mer_base_static<U>::k();
    if(start > k_ - k)
      return res;
    const size_t len = 8 * std::min(sizeof(U), sizeof(base_type));

    const unsigned int bstart = 2 * (k_ - start - k);
    for(unsigned int i = 0; (i < 2 * k) && (bstart + i < 2 * k_); i += len)
      res.set_bits(i, len, this->get_bits(bstart + i, len));

    return res;
  }

  unsigned int k() const { return k_; }
  unsigned int size() const { return k_; }
  static unsigned int k(unsigned int k) { return k; }

private:
  unsigned int k_;
};

typedef dna_fragment_base<uint64_t> dna_fragment;

#endif /* __DNA_FRAGMENT_HPP__ */
