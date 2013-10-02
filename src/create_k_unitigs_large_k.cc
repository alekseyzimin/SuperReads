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

#include <sys/types.h>
#include <signal.h>
#include <assert.h>

#include <iostream>
#include <memory>
#include <algorithm>
#include <set>

#include <src/create_k_unitigs_common.hpp>

#include <jellyfish/mapped_file.hpp>
#include <jellyfish/thread_exec.hpp>
#include <gzip_stream.hpp>
#include <src/bloom_counter2.hpp>
#include <src/MurmurHash3.h>
#include <src/read_parser.hpp>
#include <jellyfish/mer_dna.hpp>
#include <src/mer_stream.hpp>
#include <src/bloom_filter.hpp>
#include <jflib/multiplexed_io.hpp>
#include <jellyfish/atomic_field.hpp>
#include <src/create_k_unitigs_large_k_cmdline.hpp>

// GLOBAL: command line switches
cmdline_parse args;

using jellyfish::thread_exec;
using jellyfish::mer_dna;
struct mer_dna_hash {
  void operator()(const mer_dna& m, uint64_t *hashes) const {
    MurmurHash3_T_128(m, (m.k() / 4) + (m.k() % 4 != 0), 0, hashes);
  }
};
typedef bloom_counter2<mer_dna, mer_dna_hash> mer_bloom_counter2;


/* Read k-mers and store them in a map. The map_type must have the
   operator[]. The content returned must have the prefix ++. All this
   needs to be multi-thread safe.
 */
template<typename map_type, typename parser_type>
class populate_mer_set : public thread_exec {
  int          mer_len_;
  parser_type& parser_;
  map_type&    set_;

public:
  populate_mer_set(int mer_len, map_type& set, parser_type& parser) :
    mer_len_(mer_len), parser_(parser), set_(set)
  { }

  void start(int thid) {
    mer_stream<mer_dna, parser_type> stream(mer_len_, parser_);

    for( ; stream; ++stream)
      ++set_[stream.canonical()];
  }
};
typedef populate_mer_set<mer_bloom_counter2, read_parser> mer_populate;


typedef bloom_filter<mer_dna, mer_dna_hash, mt_access<unsigned int>> bloom_filter_type;
typedef create_k_unitig<mer_bloom_counter2, bloom_filter_type, bloom_filter_type, read_parser, cmdline_parse> unitiger_type;

std::ostream* open_output() {
  if(!args.output_given)
    return new std::ostream(std::cout.rdbuf());
  if(args.gzip_flag)
    return new gzipstream(args.output_arg);
  return new std::ofstream(args.output_arg);
}

int main(int argc, char *argv[])
{
  args.parse(argc, argv);

  mer_dna::k(args.mer_arg);
  if(!args.min_len_given)
    args.min_len_arg = args.mer_arg + 1;

  // Populate Bloom filter with k-mers
  std::auto_ptr<mer_bloom_counter2> kmers;
  std::auto_ptr<jellyfish::mapped_file> dbf;

  {
    if(args.load_given) {
      dbf.reset(new jellyfish::mapped_file(args.load_arg));
      uint64_t* base = (uint64_t*)dbf->base();
      kmers.reset(new mer_bloom_counter2(base[0], base[1], (unsigned char*)(base + 2)));
    } else {
      kmers.reset(new mer_bloom_counter2(args.false_positive_arg, args.nb_mers_arg));
      read_parser parser(args.input_arg.begin(), args.input_arg.end(),
                         args.threads_arg);
      mer_populate populate(args.mer_arg, *kmers, parser);
      populate.exec_join(args.threads_arg);
    }
  }

  if(args.save_given) {
    std::ofstream save_file(args.save_arg);
    if(!save_file) {
      std::cerr << "Can't open file '" << args.save_arg << "'" << std::endl;
      exit(EXIT_FAILURE);
    }
    uint64_t x;
    x = kmers->m();
    save_file.write((char*)&x, sizeof(x));
    x = kmers->k();
    save_file.write((char*)&x, sizeof(x));
    kmers->write_bits(save_file);
  }

  {
    std::auto_ptr<std::ostream> output_ostream(open_output());
    bloom_filter_type used(args.false_positive_arg, args.nb_mers_arg);
    read_parser parser(args.input_arg.begin(), args.input_arg.end(),
                       args.threads_arg);
    unitiger_type unitiger(*kmers, used, args.threads_arg, parser,
                           *output_ostream, args);
    unitiger.exec_join(args.threads_arg);
  }

  return EXIT_SUCCESS;
}
