#include <gtest/gtest.h>
#include <unittests/misc.hpp>
#include <src/read_parser.hpp>
#include <tmpstream.hpp>

static const int nb_sequences = 100000;
static const int nb_bases_dev = 100;
static const int nb_threads   = 5;

class ReadParserRandomSeq : public ::testing::Test {
protected:
  static void SetUpTestCase() {
    total_bases = 0;
    char bases[4] = { 'A', 'C', 'G', 'T' };
    for(int i = 0; i < nb_sequences; ++i) {
      random_fasta << ">" << i << "\n";
      int nb_bases = 50 + (random() % nb_bases_dev);
      for(int j = 0; j < nb_bases; ++j) {
        long r = random();
        random_fasta << bases[random() % 4];
        if((r >> 2) % 100 == 0)
          random_fasta << "\n";
      }
      random_fasta << "\n";
      total_bases += nb_bases;
    }
  }

  std::istream fasta_stream;
  int actual_bases;
public:
  ReadParserRandomSeq() :
    fasta_stream(random_fasta.rdbuf()), actual_bases(total_bases) { }

  virtual void SetUp() {
    fasta_stream.seekg(0);
    fasta_stream.clear();
  }


  static tmpstream random_fasta;
  static int total_bases;
};

tmpstream ReadParserRandomSeq::random_fasta;
int ReadParserRandomSeq::total_bases = 0;

struct sum_line_lengths_data {
  read_parser              parser;
  jflib::atomic_field<int> sum;
  sum_line_lengths_data(std::istream& is) : parser(is), sum(0) { }
};
void* sum_line_lengths(void* d) {
  auto data = (sum_line_lengths_data*)d;
  read_parser::stream read_stream(data->parser);
  int bases = 0;

  for( ; read_stream; ++read_stream)
    bases += strlen(read_stream->sequence);

  data->sum += bases;
  return 0;
}

TEST_F(ReadParserRandomSeq, Lengths) {
  sum_line_lengths_data data(fasta_stream);
  pdo(nb_threads, sum_line_lengths, (void*)&data);

  EXPECT_EQ(actual_bases, jflib::a_load(data.sum));
}
