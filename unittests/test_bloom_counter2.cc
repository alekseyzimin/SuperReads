#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <gtest/gtest.h>
#include <charb.hpp>
#include <src/bloom_counter2.hpp>

typedef std::pair<std::string,unsigned int> elt;
TEST(BloomCounter, FalsePositive) {
  static const size_t nb_inserts = 2048;
  static const double error_rate = 0.01;
  static const size_t nb_strings = 1024;
  static const size_t str_len    = 100;
  std::vector<elt>    counts;
  charb str(str_len);
  str[str_len] = '\0';

  for(size_t i = 0; i < nb_strings; ++i) {
    for(size_t j = 0; j < str_len; ++j)
      str[j] = 'A' + (random() % 26);
    counts.push_back(std::make_pair((char*)str, (unsigned int)0));
  }

  bloom_counter2<const char *> bc(error_rate, nb_inserts);
  for(size_t i = 0; i < nb_inserts; ++i) {
    std::vector<elt>::reference ref = counts[random() % nb_strings];
    ++ref.second;
    bc.insert(ref.first.c_str());
  }

  size_t nb_errors = 0;

  // Check known strings
  for(size_t i = 0; i < nb_strings; ++i) {
    std::vector<elt>::reference ref = counts[i];
    unsigned int expected = std::min(ref.second, (unsigned int)2);
    unsigned int actual = bc.check(ref.first.c_str());
    EXPECT_LE(expected, actual);
    if(expected != actual)
      ++nb_errors;
  }
  EXPECT_TRUE(error_rate * nb_strings > nb_errors);

  nb_errors = 0;
  // Check unknown strings
  for(size_t i = 0; i < nb_inserts; ++i) {
    for(size_t j = 0; j < str_len; ++j)
      str[j] = '0' + (random() % 10);
    if(bc.check(str) > 0)
      ++nb_errors;
  }
  EXPECT_TRUE(2 * error_rate * nb_inserts > nb_errors);
}
