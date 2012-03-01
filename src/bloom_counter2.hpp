#ifndef __BLOOM_COUNTER2_HPP__
#define __BLOOM_COUNTER2_HPP__

#include <math.h>
#include <limits.h>
#include <reallocators.hpp>
#include <jflib/atomic_field.hpp>
#include <jflib/compare_and_swap.hpp>
#include <src/bloom_hash.hpp>

#define LOG2    0.6931471805599453
#define LOG2_SQ 0.4804530139182014
/* Bloom counter with 3 values: 0, 1 or 2. It is thread safe and lock free.
 */
template<typename Key, typename HashPair = hash_pair<Key>, typename R=remaper<unsigned char> >
class bloom_counter2 {
  const size_t          m_;
  const unsigned long   k_;
  R                     realloc;
  unsigned char * const data_;
  HashPair              hash_fns_;

public:
  typedef Key key_type;

  bloom_counter2(double fp, size_t n) :
    m_(n * (size_t)lrint(-log(fp) / LOG2_SQ)),
    k_(lrint(-log(fp) / LOG2)),
    realloc(),
    data_(realloc(0, 0, m_ / 5 + (m_ % 5 != 0))),
    hash_fns_() { }

  bloom_counter2(size_t m, unsigned long k) :
    m_(m), k_(k), realloc(), 
    data_(realloc(0, 0, m_ / 5 + (m_ % 5 != 0))),
    hash_fns_() { }

  // Number of hash functions
  unsigned long k() const { return k_; }
  // Size of bit vector
  size_t m() const { return m_ ; }

  // Insert key k. Returns previous value of k
  unsigned int insert(const Key &k) {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return insert(hashes);
  }

  // Insert key with given hashes
  unsigned int insert(const uint64_t* hashes) {
    unsigned char res = 2;
    for(unsigned long i = 0; i < k_; ++i) {
      size_t        p    = (hashes[0] % m_ + i * (hashes[1] % m_)) % m_;
      const size_t  off  = p / 5;
      const size_t  boff = p % 5;
      unsigned char v    = jflib::a_load(&data_[off]);

      while(true) {
        unsigned char w = v;
        switch(boff) {
        case 0:          break;
        case 1: w /= 3;  break;
        case 2: w /= 9;  break;
        case 3: w /= 27; break;
        case 4: w /= 81; break;
        }
        w = w % 3;
        if(w == 2) break;
        unsigned char nv = v;
        switch(boff) {
        case 0: nv += 1;  break;
        case 1: nv += 3;  break;
        case 2: nv += 9;  break;
        case 3: nv += 27; break;
        case 4: nv += 81; break;
        }
        if(jflib::cas(&data_[off], v, nv, &v)) {
          if(w < res)
            res = w;
          break;
        }
      }
    }
    return res;
  }

  unsigned int check(const Key &k) {
    uint64_t hashes[2];
    hash_fns_(k, hashes);
    return check(hashes);
  }

  unsigned int check(uint64_t *hashes) const {
    unsigned char res = 2;
    for(unsigned int i = 0; i < k_; ++i) {
      size_t        p    = (hashes[0] % m_ + i * (hashes[1] % m_)) % m_;
      const size_t  off  = p / 5;
      const size_t  boff = p % 5;
      unsigned char w    = jflib::a_load(&data_[off]);

      switch(boff) {
      case 0:          break;
      case 1: w /= 3;  break;
      case 2: w /= 9;  break;
      case 3: w /= 27; break;
      case 4: w /= 81; break;
      }
      w = w % 3;
      if(w < res)
        res = w;
    }
    return res;
  }

  // Limited std::map interface compatibility
  class element_proxy {
    bloom_counter2& bc_;
    const Key&      k_;

  public:
    element_proxy(bloom_counter2& bc, const Key& k) : bc_(bc), k_(k) { }

    unsigned int operator++() { 
      unsigned res = bc_.insert(k_);
      return res == 0 ? 1 : 2;
    }

    unsigned int operator++(int) { return bc_.insert(k_); }
    unsigned int operator*() { return bc_.check(k_); }
    operator unsigned int() { return bc_.check(k_); }
  };
  element_proxy operator[](const Key& k) { return element_proxy(*this, k); }
};

#endif // __BLOOM_COUNTER2_HPP__
