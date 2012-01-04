#include <stdlib.h>
#include <gtest/gtest.h>
#include <src/bloom_filter.hpp>

TEST(BloomFilter, FalsePositive) {
  // Enter strings with digits only. Check presence of strings with
  // letters only.
  const unsigned long nb_entry = 1024;
  double error_rate = 0.01;
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
  ASSERT_TRUE(((double)false_positive / (double)nb_entry) < error_rate);
  
  false_positive = 0;
  for(unsigned long i = 0; i < nb_entry; ++i) {
    // Create random imput string
    for(int j = 0; j < INPUTLEN; ++j)
      input[j] = 'A' + (random() % 26);
    if(f.is_member(input))
      ++false_positive;
  }
  ASSERT_TRUE(((double)false_positive / (double)nb_entry) < 2.0 * error_rate);
}
