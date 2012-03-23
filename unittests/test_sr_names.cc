#include <limits>
#include <vector>
#include <gtest/gtest.h>
#include <src/sr_names.hpp>
#include <charb.hpp>

TEST(SR_name, encode_decode) {
  long sr_name_len = random() % 10;
  charb sr_original;
  charb one_entry;
  sprintf(sr_original, "%ld%c", random(), random() % 2 ? 'R' : 'F');
  for(long i = 1; i < sr_name_len; ++i) {
    sprintf(one_entry, "_%ld_%ld%c", random() % 100, random(), random() %2 ? 'R' : 'F');
    strcat(sr_original, one_entry);
  }

  sr_name sr(sr_original);
  charb res;
  sr.to_str(res);
  EXPECT_STREQ(sr_original, res);
}
