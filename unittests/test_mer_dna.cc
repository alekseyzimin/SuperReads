/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdio.h>
#include <gtest/gtest.h>
#include <misc.hpp>
#include <src/mer_dna.hpp>

namespace {
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
  const std::string* test_mers[4] = {
    &short_mer, &mid_mer, &long_mer, &huge_mer
  };

  TEST(KMerSimple, InitSize64) {
    struct sinfo {
      unsigned int k;
      unsigned int nb_words;
      unsigned int nb_msb;
      uint64_t     msw;
      unsigned int lshift;
    };
    sinfo ary[5] = {
      { 5, 1, 10, (uint64_t)0x3ff, 8 }, { 32, 1, 64, (uint64_t)-1, 62 }, 
      { 33, 2, 2, (uint64_t)0x3, 0 }, { 64, 2, 64, (uint64_t)-1, 62 },
      { 65, 3, 2, (uint64_t)0x3, 0 }
    };
    typedef mer_dna_ns::mer_base_dynamic<uint64_t> mer64;
    for(size_t i = 0; i < sizeof(ary) / sizeof(sinfo); ++i) {
      mer64 m(ary[i].k);
      EXPECT_EQ(ary[i].k, m.k());
      EXPECT_EQ(ary[i].nb_words, m.nb_words());
      EXPECT_EQ(ary[i].nb_msb, m.nb_msb());
      EXPECT_EQ(ary[i].lshift, m.lshift());
    }
  }

#ifdef HAVE_INT128
  TEST(KMerSimple, InitSize128) {
    struct sinfo {
      unsigned int      k;
      unsigned int      nb_words;
      unsigned int      nb_msb;
      unsigned __int128 msw;
      unsigned int      lshift;
    };
    sinfo ary[5] = {
      { 5, 1, 10, (unsigned __int128)0x3ff, 8 },
      { 32, 1, 64, (unsigned __int128)0xffffffffffffffffUL, 62 }, 
      { 33, 1, 66, (unsigned __int128)0xffffffffffffffffUL | ((unsigned __int128)0x3 << 64), 64 },
      { 64, 1, 128, (unsigned __int128)-1, 126 },
      { 65, 2, 2, (unsigned __int128)0x3, 0 }
    };
    typedef mer_dna_ns::mer_base_dynamic<unsigned __int128> mer128;
    for(size_t i = 0; i < sizeof(ary) / sizeof(sinfo); ++i) {
      mer128 m(ary[i].k);
      EXPECT_EQ(ary[i].k, m.k());
      EXPECT_EQ(ary[i].nb_words, m.nb_words());
      EXPECT_EQ(ary[i].nb_msb, m.nb_msb());
      EXPECT_EQ(ary[i].lshift, m.lshift());
    }
  }
#endif

  TEST(KMerSimple, Codes) {
    EXPECT_EQ(mer_dna::CODE_A, mer_dna::code('A'));
    EXPECT_EQ(mer_dna::CODE_A, mer_dna::code('a'));
    EXPECT_EQ(mer_dna::CODE_C, mer_dna::code('C'));
    EXPECT_EQ(mer_dna::CODE_C, mer_dna::code('c'));
    EXPECT_EQ(mer_dna::CODE_G, mer_dna::code('G'));
    EXPECT_EQ(mer_dna::CODE_G, mer_dna::code('g'));
    EXPECT_EQ(mer_dna::CODE_T, mer_dna::code('T'));
    EXPECT_EQ(mer_dna::CODE_T, mer_dna::code('t'));
    EXPECT_FALSE(mer_dna::not_dna(mer_dna::CODE_A));
    EXPECT_FALSE(mer_dna::not_dna(mer_dna::CODE_C));
    EXPECT_FALSE(mer_dna::not_dna(mer_dna::CODE_G));
    EXPECT_FALSE(mer_dna::not_dna(mer_dna::CODE_T));

    for(int c = 0; c < 256; ++c) {
      switch((char)c) {
      case 'A': case 'a':
      case 'C': case 'c':
      case 'G': case 'g':
      case 'T': case 't':
        EXPECT_FALSE(mer_dna::not_dna(mer_dna::code(c)));
        break;
      default:
        EXPECT_TRUE(mer_dna::not_dna(mer_dna::code(c)));
        break;
      }
    }
  }


  // Value Type Container class
  template <typename T, int N>
  class VTC {
  public:
    typedef T Type;
    static const int test_id = N;
  };
  template <typename T, int N>
  const int VTC<T, N>::test_id;

  template<typename VT>
  class KMer : public ::testing::Test {
  public:
    typedef typename VT::Type Type;
    void SetUp() {
      Type::k(GetParam().size());
    }
    const std::string& GetParam() const {
      return *test_mers[VT::test_id];
    }
  };
  typedef ::testing::Types<VTC<mer_dna_ns::mer_base_dynamic<uint64_t>, 0>,
                           VTC<mer_dna_ns::mer_base_dynamic<uint64_t>, 1>,
                           VTC<mer_dna_ns::mer_base_dynamic<uint64_t>, 2>,
                           VTC<mer_dna_ns::mer_base_dynamic<uint64_t>, 3>,
#ifdef HAVE_INT128
                           VTC<mer_dna_ns::mer_base_dynamic<unsigned __int128>, 0>,
                           VTC<mer_dna_ns::mer_base_dynamic<unsigned __int128>, 1>,
                           VTC<mer_dna_ns::mer_base_dynamic<unsigned __int128>, 2>,
                           VTC<mer_dna_ns::mer_base_dynamic<unsigned __int128>, 3>,
#endif
                           VTC<mer_dna_ns::mer_base_static<uint32_t>, 3>,
                           VTC<mer_dna_ns::mer_base_static<uint64_t>, 0>,
                           VTC<mer_dna_ns::mer_base_static<uint64_t>, 1>,
                           VTC<mer_dna_ns::mer_base_static<uint64_t>, 2>,
                           VTC<mer_dna_ns::mer_base_static<uint64_t>, 3>,
#ifdef HAVE_INT128
                           VTC<mer_dna_ns::mer_base_static<unsigned __int128>, 0>,
                           VTC<mer_dna_ns::mer_base_static<unsigned __int128>, 1>,
                           VTC<mer_dna_ns::mer_base_static<unsigned __int128>, 2>,
                           VTC<mer_dna_ns::mer_base_static<unsigned __int128>, 3>,
#endif
                           VTC<mer_dna_ns::mer_base_static<uint32_t>, 3>
                           > KMerTypes;
  TYPED_TEST_CASE(KMer, KMerTypes);

  TYPED_TEST(KMer, InitFromStr) {
    typename TypeParam::Type m(this->GetParam());
    EXPECT_EQ(this->GetParam().size(), m.k());
    EXPECT_EQ(this->GetParam(), m.to_str());
  }

  TYPED_TEST(KMer, ShiftLeft) {
    typename TypeParam::Type m(this->GetParam().size());
    m.polyA();
    int inserted = 0;
    for(auto it = this->GetParam().begin(); it != this->GetParam().end(); ++it, ++inserted) {
      m.shift_left(*it);
    
      int check = inserted;
      for(auto cit = this->GetParam().begin(); check >= 0; ++cit, --check)
        EXPECT_EQ(*cit, (char)m.base(check));
    }
    EXPECT_EQ(this->GetParam(), m.to_str());
  }

  TYPED_TEST(KMer, ShiftRight) {
    typename TypeParam::Type m(this->GetParam().size());
    m.polyA();
    int inserted = 0;
    for(auto it = this->GetParam().rbegin(); it != this->GetParam().rend(); ++it, ++inserted) {
      m.shift_right(*it);

      int check = inserted;
      for(auto cit = this->GetParam().rbegin(); check >= 0; ++cit, --check)
        EXPECT_EQ(*cit, (char)m.base(m.k() - 1 - check));
    }
    EXPECT_EQ(this->GetParam(), m.to_str());
  }

  TYPED_TEST(KMer, Equality) {
    typename TypeParam::Type m1(this->GetParam());
    typename TypeParam::Type m2(this->GetParam().size());

    m2.polyA();

    char str[this->GetParam().size() + 1];
    memset(str, 'A', this->GetParam().size());
    str[this->GetParam().size()] = '\0';
    EXPECT_STREQ(str, m2.to_str().c_str());
  
    int i = 1;
    for(auto it = this->GetParam().begin(); it < this->GetParam().end(); ++it, ++i) {
      sprintf(str + this->GetParam().size() - i, "%.*s", i, this->GetParam().c_str());
      typename TypeParam::Type m(str);
      EXPECT_STREQ(str, m.to_str().c_str());
      m2.shift_left(*it);
      EXPECT_EQ(m, m2);
    }
    EXPECT_TRUE(m1 == m2);
    EXPECT_FALSE(m1 != m2);
    EXPECT_EQ(m1.to_str(), m2.to_str());

    // typename TypeParam::Type m3(this->GetParam());
    // m3[0] = 0;
    // EXPECT_FALSE(m1 == m3);

    // typename TypeParam::Type m4(this->GetParam().size() + 1);
    // EXPECT_FALSE(m1 == m4);
    // typename TypeParam::Type m5(this->GetParam().size() - 1);
    // EXPECT_FALSE(m1 == m5);
  
  }

  TYPED_TEST(KMer, Copy) {
    typename TypeParam::Type m1(this->GetParam());
    typename TypeParam::Type m2(m1);
    typename TypeParam::Type m3(this->GetParam().size());
    m3 = m1;

    EXPECT_TRUE(m1 == m2);
    EXPECT_TRUE(m2 == m3);
    EXPECT_TRUE(m3 == m1);
    m1.shift_left('A');
    EXPECT_TRUE(!(m1 == m2));
    EXPECT_TRUE(!(m1 == m3));
  }

  TYPED_TEST(KMer, OperatorShift) {
    typename TypeParam::Type m(this->GetParam());
    std::ostringstream os;
    os << m;
    EXPECT_EQ(this->GetParam(), os.str());
  }

  TYPED_TEST(KMer, GetBits) {
    typename TypeParam::Type m(this->GetParam());
    for(unsigned int i = 0; i < 20; ++i) {
      long int start   = random() % (this->GetParam().size() - 1);
      long int max_len = 
        std::min(this->GetParam().size() - start, 8 * sizeof(typename TypeParam::Type::base_type));
      long int len     = (random() % (max_len - 1)) + 1;
  
      // Get bits by right-shifting
      typename TypeParam::Type cm(m);
      for(unsigned int j = 1; j < start; j += 2)
        cm.shift_right(0); // Shift by 2 bits
      typename TypeParam::Type::base_type y = cm.word(0);
      if(start & 0x1)
        y >>= 1;
      y &= ((typename TypeParam::Type::base_type)1 << len) - 1;

      EXPECT_EQ(y, m.get_bits(start, len));
    }
  }
  TYPED_TEST(KMer, GetBases) {
    typename TypeParam::Type m(this->GetParam());

    for(auto it = this->GetParam().rbegin(); it != this->GetParam().rend(); ++it)
      EXPECT_EQ(*it, (char)m.base(it - this->GetParam().rbegin()));

    const char bases[4] = { 'A', 'C', 'G', 'T' };
    for(auto it = bases; it != bases + 4; ++it) {
      typename TypeParam::Type n(m);
      for(size_t j = 0; j < this->GetParam().size(); ++j)
        n.base(j) = *it;
      typename TypeParam::Type m_expected(std::string(this->GetParam().size(), *it));
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

  TYPED_TEST(KMer, ReverseComplement) {
    typename TypeParam::Type m(this->GetParam());
    std::string rc(this->GetParam().size(), 'A');
    for(size_t i = 0; i < rc.size(); ++i)
      rc[i] = rc_base(this->GetParam()[this->GetParam().size() - 1 - i]);
    EXPECT_EQ(this->GetParam().size(), m.k());
    m.reverse_complement();
    EXPECT_EQ(rc, m.to_str());
    typename TypeParam::Type rm(rc);
    EXPECT_EQ(rm, m);
    EXPECT_EQ(m, m.get_reverse_complement().get_reverse_complement());
  }

  TYPED_TEST(KMer, Canonical) {
    typename TypeParam::Type m(this->GetParam());
    typename TypeParam::Type canonical = m.get_canonical();

    EXPECT_FALSE(m < canonical);
    EXPECT_TRUE(canonical == m || canonical == m.get_reverse_complement());
    m.canonicalize();
    EXPECT_EQ(canonical, m.get_canonical());
  }
}
