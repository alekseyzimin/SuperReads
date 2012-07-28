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
      EXPECT_EQ(r, d.remainder(a));
      uint64_t qd, rd;
      d.division(a, qd, rd);
      EXPECT_EQ(q, qd);
      EXPECT_EQ(r, rd);
  }

  TEST(Divisor, Random) {
    for(int i = 0; i < 100000; ++i) {
      uint64_t a = large_random();
      uint64_t b = large_random() & (((uint64_t)1 << 62) - 1);
      test_div(a, b);
    }
  }

  TEST(Divisor, Specific) {
    struct pair_test {
      uint64_t dividend, divisor;
    };
    pair_test tests[2] = {
      { 11529215046068486072UL, 20480UL },
      { 16563367944624636593UL, 20480UL }
    };
    for(size_t i = 0; i < sizeof(tests) / sizeof(pair_test); ++i) {
      std::ostringstream scope;
      scope << tests[i].dividend << " / " << tests[i].divisor;
      SCOPED_TRACE(scope.str().c_str());
      test_div(tests[i].dividend, tests[i].divisor);
    }
  }

  TEST(Divisor, Print) {
    const std::string res("d:2947104,p:22,m:7806571878267802146");
    std::ostringstream resos;
    divisor64 d(2947104);

    resos << d;

    EXPECT_EQ(res, resos.str());
  }
}
