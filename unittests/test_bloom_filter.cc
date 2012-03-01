#include <stdlib.h>
#include <gtest/gtest.h>
#include <src/bloom_filter.hpp>
#include <src/bloom_hash.hpp>

// XXX: this does not really check the muti-thread safe behavior

static const unsigned long nb_entry   = 1024;
static const double        error_rate = 0.01;

template<typename T>
class BloomFilter : public ::testing::Test {
public:
  BloomFilter() : f(error_rate, nb_entry) { }
  T f;
};
typedef ::testing::Types<bloom_filter<const char*>,
                         bloom_filter<const char*, hash_pair<const char*>, mt_access<unsigned int>>> bloom_filter_types;
TYPED_TEST_CASE(BloomFilter, bloom_filter_types);

TYPED_TEST(BloomFilter, FalsePositive) {
  // Enter strings with digits only. Check presence of strings with
  // letters only.
  bloom_filter<const char*> f(error_rate, nb_entry);

  // The tests comparing to the error rate can fail (it is randomized
  // after all). Hopefully it is rare enough to not happen at all.

#define INPUTLEN 100
  char input[INPUTLEN + 1];
  input[INPUTLEN] = '\0';
  size_t false_positive = 0;
  for(unsigned long i = 0; i < nb_entry; ++i) {
    // Create random imput string
    for(int j = 0; j < INPUTLEN; ++j)
      input[j] = '0' + (random() % 10);
    bool in = f.insert(input);
    if(in)
      ++false_positive;
    ASSERT_TRUE(f.is_member(input));
  }
  ASSERT_GT(error_rate, ((double)false_positive / (double)nb_entry));
  
  false_positive = 0;
  for(unsigned long i = 0; i < nb_entry; ++i) {
    // Create random imput string
    for(int j = 0; j < INPUTLEN; ++j)
      input[j] = 'A' + (random() % 26);
    if(f.is_member(input))
      ++false_positive;
  }
  ASSERT_GT(2.0 * error_rate, ((double)false_positive / (double)nb_entry));
}
