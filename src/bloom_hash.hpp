#ifndef __BLOOM_HASH_HPP__
#define __BLOOM_HASH_HPP__

#include <src/MurmurHash3.h>

/* Hash pairs
 */
template<typename Key>
class hash_pair { };

template <>
class hash_pair<const char*> {
public:
  void operator()(const char* const key, uint64_t *hashes) const {
    MurmurHash3_x64_128(key, strlen(key), 0, hashes);
  }
};

#endif // __BLOOM_HASH_HPP__
