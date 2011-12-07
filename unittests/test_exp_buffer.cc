#include <gtest/gtest.h>
#include <exp_buffer.hpp>
#include <algorithm>

typedef ExpandingBuffer<int, remaper<int> > int_buf;
typedef ExpandingBuffer<int, remaper_init<int, -1> > init_buf;

TEST(Remaper, Init) {
  int_buf b;

  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)0, b.size());

  b[5] = 5;
  EXPECT_EQ((size_t)6, b.capacity());
  EXPECT_EQ((size_t)6, b.size());
  for(int i = 0; i < 5; ++i)
    EXPECT_EQ(0, b[i]);
  EXPECT_EQ(5, b[5]);
  b[3] = 3;
  EXPECT_EQ((size_t)6, b.capacity());
  EXPECT_EQ((size_t)6, b.size());
  EXPECT_EQ(3, b[3]);

  b[6] = 6;
  EXPECT_EQ((size_t)12, b.capacity());
  EXPECT_EQ((size_t)7, b.size());

  b[5000] = 5000;
  EXPECT_EQ((size_t)5001, b.size());
  EXPECT_EQ(5000, b.back());
  for(int_buf::iterator it = b.begin(); it != b.end(); ++it)
    EXPECT_TRUE(*it == 0 || *it == (it - b.begin()));
}

TEST(Remaper, Default) {
  init_buf b;
  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)0, b.size());
  b[10] = 5;
  EXPECT_EQ(b.size() - 1, (size_t)std::count(b.begin(), b.end(), -1));
}

TEST(Remaper, Swap) {
  int_buf b(10);

  for(size_t i = 0; i < b.capacity(); ++i)
    b[i] = 2 * i;
  
  EXPECT_EQ((size_t)10, b.size());
  
  int_buf bs;
  b.swap(bs);
  EXPECT_EQ((size_t)0, b.capacity());
  EXPECT_EQ((size_t)10, bs.capacity());
  for(size_t i = 0; i < bs.capacity(); ++i)
    EXPECT_EQ((int)(2 * i), bs[i]);

  std::swap(b, bs);
  
  EXPECT_EQ((size_t)10, b.capacity());
  EXPECT_EQ((size_t)0, bs.capacity());
  for(size_t i = 0; i < b.capacity(); ++i)
    EXPECT_EQ((int)(2 * i), b[i]);

  ExpBuffer<int> b1, b2(5);
  std::swap(b1, b2);
  EXPECT_EQ((size_t)5, b1.capacity());
  EXPECT_EQ((size_t)0, b2.capacity());
}
