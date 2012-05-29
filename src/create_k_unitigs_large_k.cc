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


#include <thread_exec.hpp>
#include <src/bloom_counter2.hpp>
#include <src/MurmurHash3.h>
#include <src/read_parser.hpp>
#include <src/mer_dna.hpp>
#include <src/create_k_unitigs_large_k_cmdline.hpp>

struct mer_dna_hash {
  void operator()(const mer_dna& m, uint64_t *hashes) const {
    MurmurHash3_T_128(m, (m.k() / 4) + (m.k() % 4 != 0), 0, hashes);
  }
};
typedef bloom_counter2<mer_dna, mer_dna_hash> mer_bloom_counter2;

template<typename set_type>
class populate_mer_set : public thread_exec {
  int         mer_len_;
  read_parser parser_;
  set_type&   set_;

public:
  populate_mer_set(int mer_len, set_type& set,
                   int threads, const char* sequence_file) :
    mer_len_(mer_len), parser_(sequence_file, threads),
    set_(set)
  { }

  void start(int thid) {
    read_parser::stream read_stream(parser_);
    mer_dna             fmer(mer_len_), rmer(mer_len_);

    for( ; read_stream; ++read_stream) {
      int bases = 0;
      const char *s = read_stream->sequence;
      const char *e = read_stream->sequence.end();
      for( ; s != e; ++s) {
        uint64_t code = mer_dna::code(*s);
        if(code != mer_dna::bad_code) {
          fmer.shift_left(code);
          rmer.shift_right(mer_dna::complement(code));
          if(++bases >= mer_len_)
            ++set_[fmer < rmer ? fmer : rmer];
        } else
          bases = 0;
      }
    }
  }
};
typedef populate_mer_set<mer_bloom_counter2> mer_populate;

//template<typename 

int main(int argc, char *argv[])
{
  create_k_unitigs_large_k args(argc, argv);
  
  // Populate Bloom filter with k-mers
  mer_bloom_counter2 kmers(args.false_positive_arg, args.nb_mers_arg);
  {
    mer_populate populate(args.mer_arg, kmers, args.threads_arg, args.input_arg);
    populate.exec_join(args.threads_arg);
  }

  

  return 0;
}
