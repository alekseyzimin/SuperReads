#include <gtest/gtest.h>
#include <src/exp_buffer.hpp>

TEST(Remaper, Init) {
  typedef ExpandingBuffer<int, remaper> int_buf;
  int_buf b;

  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)0, b.len());

  b[5] = 5;
  EXPECT_EQ((size_t)6, b.capacity());
  EXPECT_EQ((size_t)6, b.len());
  for(int i = 0; i < 5; ++i)
    EXPECT_EQ(0, b[i]);
  EXPECT_EQ(5, b[5]);
  b[3] = 3;
  EXPECT_EQ((size_t)6, b.capacity());
  EXPECT_EQ((size_t)6, b.len());
  EXPECT_EQ(3, b[3]);

  b[6] = 6;
  EXPECT_EQ((size_t)12, b.capacity());
  EXPECT_EQ((size_t)7, b.len());

  b[5000] = 5000;
  EXPECT_EQ((size_t)5001, b.len());
  EXPECT_EQ(5000, *(b.ptr() - 1));
  for(int_buf::iterator it = b.begin(); it != b.end(); ++it)
    EXPECT_TRUE(*it == 0 || *it == (it - b.begin()));
}
