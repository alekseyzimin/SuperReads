#include <gtest/gtest.h>
#include <iomanip>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_INT128
#include <int128.hpp>

namespace {
  std::string last(std::string str, int chars) {
    return str.substr(str.size() - chars);
  }
  
  class TestInt128 : public ::testing::TestWithParam<std::ios::fmtflags> { 
  protected:
    void conf(std::ostream& os) {
      os.setf(GetParam());
      if(os.flags() & std::ios::adjustfield)
        os << std::setw(10) << std::setfill('x');
    }

    void test64(int64_t i) {
      __int128 x = i;
      {
        std::ostringstream os128, os64;
        conf(os128); conf(os64);

        os128 << x;
        os64 << i;
        EXPECT_EQ(os64.str(), os128.str());
      }
      {
        std::ostringstream os128, os64;
        conf(os128); conf(os64);

        os128 << std::oct << x;
        os64 << std::oct << i;
        if(i < 0) {
          EXPECT_EQ(last(os64.str(), 21), last(os128.str(), 21));
          if(os64.str()[0] == '1')
            EXPECT_EQ("37", os128.str().substr(0, 2));
          else
            EXPECT_EQ("037", os128.str().substr(0, 3));
        } else {
          EXPECT_EQ(os64.str(), os128.str());
        }
      }      
      {
        std::ostringstream os128, os64;
        conf(os128); conf(os64);

        os128 << std::hex << x;
        os64 << std::hex << i;
        if(i < 0) {
          EXPECT_EQ(last(os64.str(), 16), last(os128.str(), 16));
          EXPECT_EQ(os64.str().substr(0, 3), os128.str().substr(0, 3));
        } else {
          EXPECT_EQ(os64.str(), os128.str());
        }
      }      
    }
  };

  TEST_P(TestInt128, PrintSigned64) {
    for(int64_t i = -1023; i <= 1023; ++i)
      test64(i);

    for(int i = 10; i < 63; ++i) {
      for(int j = -10; j <= 10; ++j)
        test64(((int64_t)1 << i) + j);
    }
  }

  INSTANTIATE_TEST_CASE_P(FormatFlags, TestInt128,
                          ::testing::Values(std::ios::boolalpha, // No effect
                                            std::ios::showbase,
                                            std::ios::uppercase,
                                            std::ios::uppercase | std::ios::showbase,
                                            std::ios::showpos,
                                            std::ios::left,
                                            std::ios::right,
                                            std::ios::internal));

  void test128(std::string r, __int128 x) {
    std::ostringstream os;
    os << x;
    EXPECT_EQ(r, os.str());
  }

  void testu128(std::string r, unsigned __int128 x) {
    std::ostringstream os;
    os << x;
    EXPECT_EQ(r, os.str());
  }

  TEST(TestInt128, PrintSigned128) {
    test128("1267650600228229401496703205376", (__int128)1 << 100);
    test128("-633825300187901677043189809151", -(((__int128)1 << 66) + ((__int128)1 << 99) - (__int128)1));
  }

  TEST(TestInt128, LimitsSigned) {
    test128("170141183460469231731687303715884105727", std::numeric_limits<__int128>::max());
    test128("-170141183460469231731687303715884105728", std::numeric_limits<__int128>::min());
    EXPECT_EQ(127, std::numeric_limits<__int128>::digits);
    EXPECT_EQ(38, std::numeric_limits<__int128>::digits10);
  }

  TEST(TestInt128, LimitsUnsigned) {
    testu128("340282366920938463463374607431768211455",
            std::numeric_limits<unsigned __int128>::max());
    testu128("0", std::numeric_limits<unsigned __int128>::min());
    EXPECT_EQ(128, std::numeric_limits<unsigned __int128>::digits);
    EXPECT_EQ(39, std::numeric_limits<unsigned __int128>::digits10);
  }

}
#endif
