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

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#include <jellyfish/jellyfish.hpp>
#include <src/dna_fragment.hpp>
#include <src/filter_overlap_file_cmdline.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jflib/multiplexed_parser.hpp>
#include <jflib/multiplexed_io.hpp>
#include <jellyfish/binary_dumper.hpp>
#include <jellyfish/large_hash_array.hpp>

using jellyfish::thread_exec;
using jellyfish::mer_dna;
using jflib::o_multiplexer;
using jflib::omstream;
typedef std::unordered_map<unsigned long, dna_fragment> fragment_map;

fragment_map read_fragments(const char* path) {
  std::ifstream fd(path);
  std::string   header, sequence;
  fragment_map  res;

  while(true) {
    if(!getline(fd, header)) break;
    if(header[0] != '>') break;
    if(!getline(fd, sequence)) break;
    size_t        comma = header.find(',');
    unsigned long id    = std::stoul(header.substr(comma + 1));
    res[id]             = dna_fragment(sequence);
  }

  return res;
}

std::streampos remaining_length(std::ifstream& is) {
  std::streampos current = is.tellg();
  if(current == -1)
    return 0;
  is.seekg(0, std::ios::end);
  if(!is.good()) {
    is.clear();
    return 0;
  }
  std::streampos end = is.tellg();
  is.seekg(current, std::ios::beg);
  return end == -1 ? 0 : end - current;
}

mer_array read_bad_mers(const char* path) {
  std::ifstream fd(path);
  if(!fd.good())
    die << "Can't open hash '" << path << "'";
  jellyfish::file_header header(fd);
  // Should check format
  mer_dna::k(header.key_len() / 2);
  binary_reader reader(fd, &header);
  mer_array ary(header.size(), // TODO: estimate real number of mers
                header.key_len(),
                0, // No value associated, it is a set
                header.max_reprobe());
  while(reader.next()) {
    ary.set(reader.key());
  }

  return ary;
}

class line_parser : public multiplexed_parser<std::string> {
  typedef multiplexed_parser<std::string> super;
  std::ifstream fd;
public:
  line_parser(int nb_threads, const char* path) :
    super(nb_threads, 100),
    fd(path)
  { }

  bool good() const { return fd.good(); }

  virtual void parser_loop() {
    while(fd) {
      elt e(elt_init());
      size_type& i = e->nb_filled;
      for(i = 0; fd && i < group_size(); ++i)
        getline(fd, e->elements[i]);
    }
  }
};

class filter : public thread_exec {
  line_parser         lines_;
  const fragment_map& fragments_;
  const mer_array&    bad_mers_;
  std::ofstream       ofd_;
  o_multiplexer       mult_;

public:
  filter(int nb_threads, const char* overlaps, const fragment_map& fragments,
         const mer_array& bad_mers, const char* output) :
    lines_(nb_threads, overlaps), fragments_(fragments), bad_mers_(bad_mers),
    ofd_(output), mult_(&ofd_, 3 * nb_threads, 4096)
  {
    if(!lines_.good())
      die << "Failed to open overlap file '" << overlaps << "'";
    if(!ofd_.good())
      die << "Failed to open output file '" << output << "'";
    lines_.start_parsing();
  }

  virtual void start(int id) {
    // hashes to store k-mers from each reads
    mer_array mers1(4096, // fragment at most 2047 long
                    mer_dna::k() * 2,
                    0, // keep a set
                    62);
    mer_array mers2(4096, mer_dna::k() * 2, 0, 62);
    omstream out(mult_);

    for(line_parser::stream stream(lines_); stream; ++stream) {
      const char* line = stream->c_str();
      char*       endptr;
      // Parse line: id1 id2 ori ahang bhang ....
      unsigned long id1 = strtol(line, &endptr, 10);
      unsigned long id2 = strtol(endptr, &endptr, 10);
      if(id1 >= id2)
        continue;
      while(*endptr && isspace(*endptr)) ++endptr;
      if(*endptr == '\0')
        continue;
      char ori = *endptr++;
      long ahang = strtol(endptr, &endptr, 10);
      long bhang = strtol(endptr, &endptr, 10);
      //      std::cerr << id1 << " " << id2 << " " << ori << " " << ahang << " " << bhang << "\n";

      // Get sequences and range of matches
      auto frag1_it = fragments_.find(id1);
      if(frag1_it == fragments_.end()) {
        std::cerr << "Unknown fragment id " << id1 << "\n";
        continue;
      }
      const dna_fragment& frag1 = frag1_it->second;

      auto frag2_it = fragments_.find(id2);
      if(frag2_it == fragments_.end()) {
        std::cerr << "Unknown fragment id " << id2 << "\n";
        continue;
      }
      dna_fragment frag2 = frag2_it->second;
      if(ori == 'I')
        frag2.reverse_complement();

      long start1 = 0, start2 = 0;
      if(ahang >= 0) {
        start1 = ahang;
      } else {
        start2 = -ahang;
      }
      long end1 = frag1.size(), end2 = frag2.size();
      if(bhang >= 0) {
        end2 -= bhang;
      } else {
        end1 += bhang;
      }
      //      std::cerr << "(" << start1 << ", " << end1 << ") ("
      //                << start2 << ", " << end2 << ") " << ori << "\n";

      // Find k-mers in common in the overlap region
      mers1.clear();
      for( ; start1 <= end1 - mer_dna::k(); ++start1) {
        mer_dna m1 = frag1.sub_mer(start1);
        //        std::cerr << "mer1 " << m1 << "\n";
        mers1.set(m1);
      }

      mers2.clear();
      for( ; start2 <= end2 - mer_dna::k(); ++start2) {
        mer_dna m2 = frag2.sub_mer(start2);
        //        std::cerr << "mer2 " << m2 << " in 1 " << mers1.has_key(m2) << "\n";
        if(mers1.has_key(m2))
          mers2.set(m2);
      }

      auto it = mers2.eager_slice(0, 1);
      while(it.next()) {
        //        std::cerr << "test " << it.key() << (bad_mers_.has_key(it.key().get_canonical()) ? " bad" : " good") << "\n";
        if(!bad_mers_.has_key(it.key().get_canonical())) {
          out << id1 << " " << id2 << " " << ori << " " << ahang << " " << bhang << endptr << "\n";
          out << jflib::endr;
          break;
        }
      }
    }
  }
};

int main(int argc, char *argv[])
{
  filter_overlap_file_cmdline args(argc, argv);

  fragment_map fragments = read_fragments(args.fragments_arg);
  mer_array bad_mers = read_bad_mers(args.kmer_arg);

  filter overlap_filter(args.threads_arg, args.overlaps_arg, fragments, bad_mers, args.output_arg);
  overlap_filter.exec_join(args.threads_arg);

  return 0;
}
