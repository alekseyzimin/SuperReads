#include <gtest/gtest.h>
#include <jellyfish/mer_dna.hpp>
#include <src/dna_fragment.hpp>
#include <misc.hpp>

namespace {
using jellyfish::mer_dna;

TEST(DNAFragment, SubMers) {
  mer_dna::k(random_bits(5) + 20);
  const size_t flength = random_bits(12) + mer_dna::k() + 100;
  std::string frag_seq;
  for(size_t i = 0; i < flength; ++i)
    frag_seq += mer_dna::rev_code(random_bits(2));
  dna_fragment fragment(frag_seq);
  EXPECT_EQ(frag_seq, fragment.to_str());

  for(size_t i = 0; i <= flength - mer_dna::k(); ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i << " flength:" << flength 
                 << " k:" << mer_dna::k());
    mer_dna m = fragment.sub_mer<uint64_t>(i);
    EXPECT_EQ(frag_seq.substr(i, mer_dna::k()), m.to_str());
  }
}
}
