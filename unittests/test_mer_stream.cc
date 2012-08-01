#include <gtest/gtest.h>
#include <src/mer_dna.hpp>
#include <src/read_parser.hpp>
#include <src/mer_stream.hpp>

namespace {
  TEST(MerStream, Count) {
    std::string fasta(">r0\n"
                      "TGCAATGGGATTATGTCCTA\n"
                      ">r1\n"
                      "GCCACACTGGTCTGGACGACACATACAACGTCCGGTCCATAAATCCGAGAA\n"
                      ">r2\n"
                      "TCGNGCAAGGCTGATGCATCTTTATTGGCGTAGGAATGGTCACCGTTCATT\n");
    std::istringstream fasta_stream(fasta);
    read_parser parser(fasta_stream, 1);
    mer_dna::k(10);
    mer_stream<mer_dna, read_parser> stream(10, parser);
    
    int count = 0;
    for( ; stream; ++stream)
      ++count;
    EXPECT_EQ(11 + 42 + 38, count);
  }

  TEST(MerStream, Mers) {
    std::string sequence("GTTCAACTGTGCTACCCCATGGGTGGCAGGCGCGGCGTGTCCCATATCGA");
    std::string fasta(">r0\n");
    fasta += sequence;
    std::istringstream fasta_stream(fasta);

    read_parser parser(fasta_stream, 1);
    mer_dna::k(10);
    mer_stream<mer_dna, read_parser> stream(10, parser);

    for(int i = 0; stream; ++stream, ++i) {
      EXPECT_EQ(sequence.substr(i, 10), stream.fmer().to_str());
      mer_dna rmer(stream.rmer());
      rmer.reverse_complement();
      EXPECT_EQ(sequence.substr(i, 10), rmer.to_str());
    }
  }
}
