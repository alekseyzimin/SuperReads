#ifndef __SKIP_LIST_COMMON_HPP__
#define __SKIP_LIST_COMMON_HPP__

#include <gcc_builtins.hpp>

// XOR RNG by George Marsalia
struct xor_random {
  typedef uint64_t rand_type;
  uint64_t x;
  // xor64
  rand_type operator()() {
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return x;
  }

  // xor128
  // uint64_t x, y, z, w;
  // rand_type operator()() {
  //   uint64_t t = x ^ (x << 5);
  //   x = y;
  //   y = z;
  //   z = w;
  //   w = (w ^ (w >> 29)) ^ (t ^ (t >> 12));
  //   return w;
  // }
  explicit xor_random() : x(88172645463325252LL) { };
  explicit xor_random(uint64_t seed, int n = 10) : x(seed) {
    for(int i = 0; i < n; ++i)
      this->operator()();
  }
};

// Return the height of a tower of pointer. Specialized for 2 and 4.
template<typename Random, int p>
struct random_height;
template<typename Random>
struct random_height<Random, 2> {
  Random rng;
  int operator()() {
    typename Random::rand_type x = rng();
    return (x == 0 ? 8*sizeof(typename Random::rand_type) : ctz(x)) + 1;
  }
  random_height(const Random& rng_ = Random()) : rng(rng_) { }
};
template<typename Random>
struct random_height<Random, 4> {
  Random rng;
  int operator()() {
    typename Random::rand_type x = rng();
    return (x == 0 ? 4 * sizeof(typename Random::rand_type) : ctz(x) >> 1) + 1;
  }
  random_height(const Random& rng_ = Random()) : rng(rng_) { }
};

// Comparator for the first element of a pair
template<typename Pair, class Compare>
struct first_comp {
  typedef Pair pair_type;
  typedef Compare comp_type;
  Compare comp;
  bool operator()(const pair_type& p1, const pair_type& p2) const {
    return comp(p1.first, p2.first);
  }
  bool operator()(const typename pair_type::first_type& p1, const pair_type& p2) const {
    return comp(p1, p2.first);
  }
  bool operator()(const pair_type& p1, const typename pair_type::first_type& p2) const {
    return comp(p1.first, p2);
  }
  first_comp(const Compare& comp_ = Compare()) : comp(comp_) { }
};


#endif /* __SKIP_LIST_COMMON_HPP__ */
