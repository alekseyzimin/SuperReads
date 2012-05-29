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


#ifndef __BLOOM_FILTER_HPP__
#define __BLOOM_FILTER_HPP__

#include <math.h>
#include <vector>
#include <algorithm>
#include <reallocators.hpp>
#include <src/bloom_hash.hpp>

/* Bloom filter using Kirsh & Mitzenmacher double hashing. I.e., only
   two hash functions are computed and the k functions have values (0
   <= i < k):

   hash_0 + i * hash_1
 */

/* Memory operator for serial access.
 */
template<typename T>
struct serial_access {
  static T fetch_and_or(T* ptr, T x) {
    T res = *ptr;
    *ptr |= x;
    return res;
  }

  static T fetch(T* ptr) {
    return *ptr;
  }
};

/* Memory operator for muti-threaded access.
 */
template<typename T>
struct mt_access {
  static T fetch_and_or(T* ptr, T x) { return __sync_ftetch_and_or(ptr, x); }
  static T fetch(T* ptr) { return *(volatile T*)ptr; }
};

#define LOG2    0.6931471805599453
#define LOG2_SQ 0.4804530139182014
template<typename Key, typename HashPair = hash_pair<Key>, typename M=serial_access<unsigned int>, typename R=remaper<unsigned int>>
class bloom_filter {
  typedef typename R::element_type  element_type;
  typedef typename R::element_type* element_pointer;
  const size_t        m_;
  const unsigned long k_;       // Number of hashes
  R                   realloc_;
  element_pointer     data_;
  const HashPair      hash_fns_;
  M                   mem_access_;
  
  typedef std::vector<bool>::reference bit;

  static const size_t elt_size = sizeof(element_type) * 8;
  static size_t nb_elements(size_t m) { return m / elt_size + (m % elt_size != 0); }

public:
  // BF with false positive rate of fp and estimated number of entries
  // of n.
  bloom_filter(double fp, size_t n) : 
    m_(n * (size_t)lrint(-log(fp) / LOG2_SQ)),
    k_(lrint(-log(fp) / LOG2)),
    realloc_(),
    data_(realloc_(0, 0, nb_elements(m_))),
    hash_fns_(), mem_access_()
  { }

  bloom_filter(size_t m, unsigned long k) :
    m_(m), k_(k), realloc_(),
    data_(realloc_(0, 0, nb_elements(m_))),
    hash_fns_(), mem_access_() { }

  // Number of hash functions
  unsigned long k() const { return k_; }
  // Size of bit vector
  size_t m() const { return m_; }

  // Insert key k
  bool insert(const Key &k) {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return insert(hashes);
  }
    
  // Insert key with given hashes
  bool insert(const uint64_t *hashes) {
    bool     present = true;
    
    for(unsigned long i = 0; i < k_; ++i) {
      size_t pos    = (hashes[0] + i * hashes[1]) % m_;
      size_t elt_i  = pos / elt_size;
      size_t bit_i  = pos % elt_size;
      present      &= mem_access_.fetch_and_or(data_ + elt_i, ((element_type)1) << bit_i);
    }

    return present;
  }

  // Compute hashes of key k
  void hash(const Key &k, uint64_t *hashes) const { hash_fns_(k, hashes); }

  // True if k is a member of the set
  bool is_member(const Key &k) const {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return is_member(hashes);
  }

  bool is_member(const uint64_t *hashes) const {
    for(unsigned long i = 0; i < k_; ++i) {
      size_t pos   = (hashes[0] + i * hashes[1]) % m_;
      size_t elt_i = pos / elt_size;
      size_t bit_i = pos % elt_size;
      if(!(mem_access_.fetch(data_ + elt_i) & ((element_type)1 << bit_i)))
        return false;
    }
    return true;
  }
};

#endif // __BLOOM_FILTER_HPP__
