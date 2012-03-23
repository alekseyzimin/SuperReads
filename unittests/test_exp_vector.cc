#include <gtest/gtest.h>
#include <exp_vector.hpp>

TEST(ExpVector, Subscript) {
  exp_vector<int> ev;

  ev[10] = 5;
  for(int i = 0; i < 10; ++i)
    EXPECT_EQ(0, ev[i]);
  EXPECT_EQ(5, ev[10]);
}
