#include <stdio.h>
#include <gtest/gtest.h>
#include <misc.hpp>
#include <src/mer_dna.hpp>

const std::string short_mer("ACGTTGGCCAA");
const std::string mid_mer("AAGACTTAGAATCAGCTAAGGTAACTTACGTAGAATATAGA");
const std::string long_mer("AAGACTTAGAATCAGCTAAGGTAACTTACGTAGAATATAG"
                           "AAGACGGCGCCTATCCCAACCTCCATTTGCTGACGCCCTC"
                           "GTAACCGTGCTGCGAGGTTACTCTATACTGA");
const std::string huge_mer("CAGGACACAATACGTCGATCCAATCCCCGACGTGAGAGTT"
                           "TAACGCTAATCTTGATCCATTACAAGTATAGATATTTCGG"
                           "GCGCCACGGGGAAACTTGCCCTATGTTCGAGTCCGCCACC"
                           "GCGGCACCAGCCTTTGTGTGACGGCCACAAAGGGTAAAAG"
                           "ATGTTGTCTCGCGCCGGCGTGCGGTCTTTACCAGATCCTT"
                           "GAGCGGCCTAGAAAGTTCGTACCAGTTCGACTGAAACAAG"
                           "ACAACGGATCGCCCACGGTCATCACACGCCCACGCACGGG"
                           "GTCGGGTTGGTATATTCAACCTCGAGTTAAACGT");

TEST(KMerSimple, InitSize) {
  struct s {
    unsigned int k;
    unsigned int nb_words;
    unsigned int nb_msb;
    uint64_t     msw;
    unsigned int lshift;
  };
  s ary[5] = {
    { 5, 1, 10, (uint64_t)0x3ff, 8 }, { 32, 1, 64, (uint64_t)-1, 62 }, 
    { 33, 2, 2, (uint64_t)0x3, 0 }, { 64, 2, 64, (uint64_t)-1, 62 },
    { 65, 3, 2, (uint64_t)0x3, 0 }
  };
  for(size_t i = 0; i < sizeof(ary) / sizeof(s); ++i) {
    mer_dna m(ary[i].k);
    EXPECT_EQ(ary[i].k, m.k());
    EXPECT_EQ(ary[i].nb_words, m.nb_words());
    EXPECT_EQ(ary[i].nb_msb, m.nb_msb());
    EXPECT_EQ(ary[i].lshift, m.lshift());
  }
}

class KMer : public ::testing::TestWithParam<std::string> {};
TEST_P(KMer, InitFromStr) {
  mer_dna m(GetParam());
  EXPECT_EQ(GetParam().size(), m.k());
  EXPECT_EQ(GetParam(), m.to_str());
}

TEST_P(KMer, ShiftLeft) {
  mer_dna m(GetParam().size());
  m.polyA(); // Avoid valgrind error with unizitialized values
  for(auto it = GetParam().begin(); it != GetParam().end(); ++it)
    m.shift_left(*it);
  EXPECT_EQ(GetParam(), m.to_str());
}

TEST_P(KMer, ShiftRight) {
  mer_dna m(GetParam().size());
  m.polyA(); // Avoid valgrind error with unizitialized values
  for(auto it = GetParam().rbegin();
      it != GetParam().rend(); ++it) {
    m.shift_right(*it);
  }
  EXPECT_EQ(GetParam(), m.to_str());
}

TEST_P(KMer, Equality) {
  mer_dna m1(GetParam());
  mer_dna m2(GetParam().size());

  m2.polyA();

  char str[GetParam().size() + 1];
  memset(str, 'A', GetParam().size());
  str[GetParam().size()] = '\0';
  EXPECT_STREQ(str, m2.to_str().c_str());
  
  int i = 1;
  for(auto it = GetParam().begin(); it < GetParam().end(); ++it, ++i) {
    sprintf(str + GetParam().size() - i, "%.*s", i, GetParam().c_str());
    mer_dna m(str);
    EXPECT_STREQ(str, m.to_str().c_str());
    m2.shift_left(*it);
    EXPECT_TRUE(m == m2);
  }
  EXPECT_TRUE(m1 == m2);
  EXPECT_FALSE(m1 != m2);
  EXPECT_EQ(m1.to_str(), m2.to_str());

  // mer_dna m3(GetParam());
  // m3[0] = 0;
  // EXPECT_FALSE(m1 == m3);

  mer_dna m4(GetParam().size() + 1);
  EXPECT_FALSE(m1 == m4);
  mer_dna m5(GetParam().size() - 1);
  EXPECT_FALSE(m1 == m5);
  
}

TEST_P(KMer, Copy) {
  mer_dna m1(GetParam());
  mer_dna m2(m1);
  mer_dna m3(GetParam().size());
  m3 = m1;

  EXPECT_TRUE(m1 == m2);
  EXPECT_TRUE(m2 == m3);
  EXPECT_TRUE(m3 == m1);
  m1.shift_left('A');
  EXPECT_TRUE(!(m1 == m2));
  EXPECT_TRUE(!(m1 == m3));
}

TEST_P(KMer, OperatorShift) {
  mer_dna m(GetParam());
  std::ostringstream os;
  os << m;
  EXPECT_EQ(GetParam(), os.str());
}

TEST_P(KMer, GetBits) {
  mer_dna m(GetParam());
  for(unsigned int i = 0; i < 20; ++i) {
    long int start   = random() % (GetParam().size() - 1);
    long int max_len = (GetParam().size() - start) > 64 ? 64 : GetParam().size() - start;
    long int len     = (random() % (max_len - 1)) + 1;
  
    // Get bits by right-shifting
    mer_dna cm(m);
    for(unsigned int j = 1; j < start; j += 2)
      cm.shift_right((uint64_t)0); // Shift by 2 bits
    uint64_t y = cm[0];
    if(start & 0x1)
      y >>= 1;
    y &= ((uint64_t)1 << len) - 1;

    EXPECT_EQ(y, m.get_bits(start, len));
  }
}
TEST_P(KMer, GetBases) {
  mer_dna m(GetParam());

  for(auto it = GetParam().rbegin(); it != GetParam().rend(); ++it)
    EXPECT_EQ(*it, (char)m.base(it - GetParam().rbegin()));

  const char bases[4] = { 'A', 'C', 'G', 'T' };
  for(auto it = bases; it != bases + 4; ++it) {
    mer_dna n(m);
    for(size_t j = 0; j < GetParam().size(); ++j)
      n.base(j) = *it;
    mer_dna m_expected(std::string(GetParam().size(), *it));
    EXPECT_EQ(m_expected, n);
  }
}

char rc_base(char c) {
  switch(c) {
  case 'A': case 'a': return 'T';
  case 'C': case 'c': return 'G';
  case 'G': case 'g': return 'C';
  case 'T': case 't': return 'A';
  }
  return 'A'; // Should never be reached
}

TEST_P(KMer, ReverseComplement) {
  mer_dna m(GetParam());
  std::string rc(GetParam().size(), 'A');
  for(size_t i = 0; i < rc.size(); ++i)
    rc[i] = rc_base(GetParam()[GetParam().size() - 1 - i]);

  m.reverse_complement();
  EXPECT_EQ(rc, m.to_str());
  mer_dna rm(rc);
  EXPECT_EQ(rm, m);
}

// TEST_P(KMer, SetBits) {
//   mer_dna m(GetParam());
//   for(unsigned int i = 0; i < 20; ++i) {
//     long int start   = random() % (GetParam().size() - 1);
//     long int max_len = (GetParam().size() - start) > 64 ? 64 : GetParam().size() - start;
//     long int len     = (random() % (max_len - 1)) + 1;
//     uint64_t bits = g_rand<uint64_t>();
//     m.set_bits(start, len, bits);
//     EXPECT_EQ(bits & (((uint64_t)1 << len) - 1), m.get_bits(start, len));
//   }
// }

INSTANTIATE_TEST_CASE_P(KMerFromStrings, KMer, ::testing::Values(short_mer, mid_mer, long_mer, huge_mer));
