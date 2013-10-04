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

#include <assert.h>

#include <iostream>
#include <memory>
#include <algorithm>
#include <set>

#include <src/MurmurHash3.h>

#include <jellyfish/mapped_file.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/atomic_field.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>

#include <src/create_k_unitigs_common.hpp>
#include <gzip_stream.hpp>
#include <jflib/multiplexed_io.hpp>
#include <src/create_k_unitigs_large_k_cmdline.hpp>

using jellyfish::thread_exec;
using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;
using jellyfish::mer_dna_bloom_filter;

typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager;
typedef jellyfish::mer_overlap_sequence_parser<stream_manager> sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;

struct mer_dna_hash {
  void operator()(const mer_dna& m, uint64_t *hashes) const {
    MurmurHash3_T_128(m, (m.k() / 4) + (m.k() % 4 != 0), 0, hashes);
  }
};


// GLOBAL: command line switches
cmdline_parse args;

// Wrapper around mer_iterator to satisfy common interface
class read_mers {
  mer_iterator stream_;

public:
  read_mers(sequence_parser& parser, int id) : stream_(parser, true) { }
  operator bool() const { return (void*)stream_ != 0; }
  const mer_dna* operator->() const { return stream_.operator->(); }
  const mer_dna& operator*() const { return stream_.operator*(); }
  read_mers& operator++() {
    ++stream_;
    return *this;
  }
};

// Wrapper around bloom filter class to have set compatible insert
// method.
class mer_bloom {
  jellyfish::bloom_filter<mer_dna, mer_dna_hash> bf_;

public:
  mer_bloom(double fp, size_t size) : bf_(fp, size) { }
  std::pair<unsigned int, bool> insert(const mer_dna& m) {
    unsigned int r = bf_.insert(m);
    return std::make_pair(r, r == 0);
  }
};

/* Read k-mers and store them in a map. The map_type must have the
   operator[]. The content returned must have the prefix ++. All this
   needs to be multi-thread safe.
 */
template<typename map_type, typename parser_type, typename stream_type>
class populate_mer_set : public thread_exec {
  int          mer_len_;
  parser_type& parser_;
  map_type&    set_;

public:
  populate_mer_set(int mer_len, map_type& set, parser_type& parser) :
    mer_len_(mer_len), parser_(parser), set_(set)
  { }

  void start(int thid) {
    for(stream_type stream(parser_, thid); stream; ++stream)
      ++set_[*stream];
  }
};
typedef populate_mer_set<jellyfish::bloom_counter2<mer_dna, mer_dna_hash>, sequence_parser, read_mers> mer_populate;

typedef jellyfish::bloom_counter2<mer_dna, mer_dna_hash> mer_counter;

typedef create_k_unitig<mer_counter, mer_bloom, mer_bloom,
                        sequence_parser, read_mers,
                        cmdline_parse> unitiger_type;

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
  std::auto_ptr<mer_counter> kmers;
  std::auto_ptr<jellyfish::mapped_file> dbf;

  {
    // if(args.load_given) {
    //   dbf.reset(new jellyfish::mapped_file(args.load_arg));
    //   uint64_t* base = (uint64_t*)dbf->base();
    //   kmers.reset(new mer_dna_bloom_counter(base[0], base[1], (unsigned char*)(base + 2)));
    // } else {
      kmers.reset(new mer_counter(args.false_positive_arg, 2 * args.nb_mers_arg));
      stream_manager manager(args.input_arg.begin(), args.input_arg.end());\
      sequence_parser parser(mer_dna::k(), manager.nb_streams(), 3 * args.threads_arg, 4096, manager);
      mer_populate populate(args.mer_arg, *kmers, parser);
      populate.exec_join(args.threads_arg);
      //    }
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
    mer_bloom used(args.false_positive_arg, 2 * args.nb_mers_arg);
    stream_manager manager(args.input_arg.begin(), args.input_arg.end());
    sequence_parser parser(mer_dna::k(), manager.nb_streams(), 3 * args.threads_arg, 4096, manager);
    unitiger_type unitiger(*kmers, used, args.threads_arg, parser,
                           *output_ostream, args);
    unitiger.exec_join(args.threads_arg);
  }

  return EXIT_SUCCESS;
}
