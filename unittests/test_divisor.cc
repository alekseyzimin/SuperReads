#include <stdlib.h>
#include <gtest/gtest.h>
#include <divisor.hpp>

namespace {
  uint64_t large_random() {
    return ((uint64_t)random() << 32) ^ (uint64_t)random();
  }

  void test_div(uint64_t a, uint64_t b) {
      divisor64 d(b);
      uint64_t q = a / b;
      uint64_t r = a % b;

      EXPECT_EQ(q, d.divide(a));
      EXPECT_EQ(q, a / d);
      EXPECT_EQ(r, d.remainder(a));
      EXPECT_EQ(r, a % d);
      uint64_t qd, rd;
      d.division(a, qd, rd);
      EXPECT_EQ(q, qd);
      EXPECT_EQ(r, rd);
  }

  TEST(Divisor, Random) {
    static int mi = 10000000;
    for(int i = 0; i < mi; ++i) {
      uint64_t a = large_random();
      uint64_t b = large_random();
      test_div(a, b);
    }
    for(int i = 0; i < mi; ++i) {
      uint64_t a = large_random();
      uint64_t b = large_random() / (i + 1);
      test_div(a, b);
    }
  }

  TEST(Divisor, Specific) {
    struct pair_test {
      uint64_t dividend, divisor;
    };
    pair_test tests[3] = {
      { 11529215046068486072UL, 20480UL },
      { 16563367944624636593UL, 20480UL },
      { 5UL, 2UL }
    };
    for(size_t i = 0; i < sizeof(tests) / sizeof(pair_test); ++i) {
      SCOPED_TRACE(::testing::Message() << tests[i].dividend << " / " << tests[i].divisor);
      test_div(tests[i].dividend, tests[i].divisor);
    }
  }

  TEST(Divisor, Print) {
#ifdef HAVE_INT128
    const std::string res("d:2947104,p:22,m:7806571878267802146");
#else
    const std::string res("d:2947104,p:0,m:0");
#endif
    std::ostringstream resos;
    divisor64 d(2947104);

    resos << d;

    EXPECT_EQ(res, resos.str());
  }
}
