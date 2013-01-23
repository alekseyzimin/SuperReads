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
#include <jellyfish/mapped_file.hpp>
#include <thread_exec.hpp>
#include <gzip_stream.hpp>
#include <src/bloom_counter2.hpp>
#include <src/MurmurHash3.h>
#include <src/read_parser.hpp>
#include <src/mer_dna.hpp>
#include <src/mer_stream.hpp>
#include <src/bloom_filter.hpp>
#include <jflib/multiplexed_io.hpp>
#include <jflib/atomic_field.hpp>
#include <src/create_k_unitigs_large_k_cmdline.hpp>

// GLOBAL: command line switches
cmdline_parse args;

typedef mer_dna_ns::mer_base_dynamic<uint64_t> mer_type;

struct mer_dna_hash {
  void operator()(const mer_type& m, uint64_t *hashes) const {
    MurmurHash3_T_128(m, (m.k() / 4) + (m.k() % 4 != 0), 0, hashes);
  }
};
typedef bloom_counter2<mer_type, mer_dna_hash> mer_bloom_counter2;


/* Read k-mers and store them in a map. The map_type must have the
   operator[]. The content returned must have the prefix ++. All this
   needs to be multi-thread safe.

   It populates two database of k-mers. The first one is for long
   k-mers. The second one is for small k-mers which are the prefix and
   suffix of the longer k-mers. For every long mer m which is seen for
   the first time, add 1 to the count of the prefix and suffix short
   mers.
 */
template<typename map_type, typename parser_type>
class populate_mer_set : public thread_exec {
  int          mer_len_;
  int          short_len_;
  parser_type& parser_;
  map_type&    set_;
  map_type&    set_short_;

public:
  populate_mer_set(int mer_len, int short_len, map_type& set, map_type& set_short, parser_type& parser) :
    mer_len_(mer_len), short_len_(short_len), parser_(parser), set_(set), set_short_(set_short)
  { }

  void start(int thid) {
    mer_stream<mer_type, parser_type> stream(mer_len_, parser_);
    mer_type m(short_len_);

    for( ; stream; ++stream) {
      unsigned int count = ++set_[stream.canonical()];
      if(count == 1) {
        m = stream.fmer();
        m.canonicalize();
        ++set_short_[m];
        m = stream.rmer();
        m.canonicalize();
        ++set_short_[m];
      }
    }
  }
};
typedef populate_mer_set<mer_bloom_counter2, read_parser> mer_populate;

// Insert a mer in a set and return true if the k-mer is new in the
// set.
template<typename set_type>
bool insert_canonical(set_type& set, const mer_type& mer) {
  return set.insert(mer.get_canonical()).second;
  // mer_type rc(mer);
  // rc.reverse_complement();
  // bool res = set.insert(rc < mer ? rc : mer).second;
  // return res;
}

/* - mer_counts_type maps k-mer to counts. Has operator[] returning
     the count.

   - used_type is a set type with operator insert.

   - end_points_type is a set type with operator insert.

   - parser_type is the type of sequence parser
 */
template<typename mer_counts_type, typename used_type, typename end_points_type, typename parser_type>
class create_k_unitig : public thread_exec {
  const mer_counts_type&        counts_; // Counts for k-mers
  used_type                     used_mers_; // Mark all k-mers whether they have been visited already
  end_points_type               end_points_; // End points of k-unitigs, to ouput only once
  int                           threads_;
  parser_type&                  parser_;
  jflib::o_multiplexer          output_multiplexer_;
  jflib::atomic_field<uint64_t> unitig_id_;

  enum direction { forward = 1, backward = -1 };
  static direction rev_direction(direction dir) { return (direction)-dir; }

public:
  create_k_unitig(const mer_counts_type& counts, used_type& used,
                  int threads, parser_type& parser, std::ostream& output) :
    counts_(counts), used_mers_(used),
    end_points_(args.false_positive_arg, args.nb_mers_arg / 10),
    threads_(threads), parser_(parser), output_multiplexer_(&output, 3 * threads, 4096),
    unitig_id_(0)
  { }

  virtual void start(int thid) {
    mer_stream<mer_type, read_parser> stream(args.mer_arg, parser_);
    jflib::omstream                  output(output_multiplexer_);
    mer_type                          current(args.mer_arg);
    mer_type                          continuation(args.mer_arg);
    mer_type                          tmp(args.mer_arg);

    for( ; stream; ++stream) {
      auto is_new = used_mers_.insert(stream.canonical());
      if(!is_new.second)
        continue;
      // Never start a unitig on low count
      if(counts_[stream.canonical()] < args.quality_threshold_arg)
        continue;
      current = *stream;

      // Grow unitig if a starting (branching) mer
      if(starting_mer(forward, current)) {
        grow_unitig(backward, current, output);
      } else if(starting_mer(backward, current)) {
        grow_unitig(forward, current, output);
      }
      // Unique continuation on both sides -> middle of k-unitig: do nothing
    }
  }

private:
  // Check all k-mers extending in one direction. If unique
  // continuation, store it in cont and return true. A unique
  // continuation may have a count less than the min quality
  // threshold, as will be reported in the count output argument.
  //
  // On the other hand, continuation with count less than the min
  // quality threshold do not create a branch compared to a high
  // quality continuation. I.e., if the threshold is 2, and the counts
  // are as follow:
  //
  // 0, 1, 0, 0 -> unique continuation, count of 1
  // 0, 1, 0, 1 -> no unique continuation
  // 2, 0, 0, 0 -> unique continuation, count of 2
  // 2, 1, 0, 0 -> unique continuation, count of 2
  // 2, 0, 3, 0 -> no unique continuation
  //
  // Otherwise return false and the value of cont is undetermined. If
  // true is returned and count is not NULL, the count of the unique
  // continuation mer is stored in the pointer.
  bool next_mer(const direction dir, const mer_type& start, mer_type& cont,
                unsigned int* count = 0) {
    int     index;
    mer_type cont_comp(start);
    cont_comp.reverse_complement();
    cont = start;

    if(dir == forward) {
      cont.shift_left(0);
      cont_comp.shift_right(0);
      index = 0;
    } else {
      cont.shift_right(0);
      cont_comp.shift_left(0);
      index = cont.k() - 1;
    }
    auto base = cont.base(index); // Point to first or last base. Correct base to change
    auto base_comp = cont_comp.base(cont.k() - 1 - index);

    int          nb_cont    = 0, nb_low_cont = 0;
    int          code       = 0, low_code = 0;
    unsigned int cont_count = 0, low_cont_count = 0;
    for(int i = 0; i < 4; ++i) {
      base      = i;
      base_comp = mer_type::complement(i);

      unsigned int cont_count_ = counts_[cont < cont_comp ? cont : cont_comp];
      if(cont_count_ >= args.quality_threshold_arg) {
        if(++nb_cont > 1)
          return false;
        code       = i;
        cont_count = cont_count_;
      } else if(cont_count_ > 0) {
        ++nb_low_cont;
        low_code       = i;
        low_cont_count = cont_count_;
      }
    }

    if(nb_cont == 1) {
      base = code;
      if(count)
        *count = cont_count;
      return true;
    } else if(nb_cont == 0 && nb_low_cont == 1) {
        base = low_code;
        if(count)
          *count = low_cont_count;
        return true;
    }
    return false;
  }

  // Return true if m is a starting mer in the given dir: a mer is a
  // starting mer if it is branching forward, or backward, or dries
  // out, maybe after some number of low count mer to skip.
  bool starting_mer(direction dir, mer_type m) {
    int     low_cont = args.cont_on_low_arg;
    mer_type tmp1(args.mer_arg), tmp2(args.mer_arg);

    while(true) {
      unsigned int count = 0;
      if(!next_mer(dir, m, tmp1, &count))
        return true;
      if(count >= args.quality_threshold_arg) {
        if(!next_mer(rev_direction(dir), tmp1, tmp2, &count))
          return true;
        break;
      }
      if(--low_cont < 0)
        return true;
      m = tmp1;
    }
    return false;
  }

  void grow_unitig(const direction dir, const mer_type& start, jflib::omstream& output) {
    bool start_new = insert_canonical(end_points_, start);
    if(!start_new)
      return;

    mer_type            mer1(start);
    mer_type            mer2(args.mer_arg);
    mer_type            mer3(args.mer_arg);
    mer_type           *current = &mer1;
    mer_type           *cont    = &mer2;
    unsigned int       count   = 0;
    unsigned int       low_run = 0;
    unsigned int       index   = dir == forward ? 0 : start.k() - 1;
    std::string        seq;
    std::set<mer_type>  set; // Set of used mers to avoid endless loop

    while(true) {
      insert_canonical(used_mers_, *current);
      if(!insert_canonical(set, *current))
        return; // loop. Don't output anything
      if(!next_mer(dir, *current, *cont, &count))
        break;
      if(!next_mer(rev_direction(dir), *cont, mer3))
        break;
      // This can happen (only) with continuation on low. It does not
      // create a branch as far as next_mer is concerned if one low
      // count and one high count, but it still a branch in this case:
      // there are two way to go through that region
      if(mer3 != *current)
        break;
      seq += (char)cont->base(index);

      if(count < args.quality_threshold_arg) {
        if(++low_run > args.cont_on_low_arg)
          break;
      } else
        low_run = 0;

      std::swap(current, cont);
    }

    // Erase trailing low quality bases if any and reset current to be
    // the actual last k-mer (with only good quality bases). Needed
    // for the test of already written k-unitigs to be accurate.
    if(low_run > 0) {
      seq.erase(seq.size() - std::min((unsigned int)seq.size(), low_run));
      if(seq.size() >= current->k()) {
        if(dir == forward) {
          *current = seq.substr(seq.size() - current->k());
        } else {
          std::string end = seq.substr(seq.size() - current->k());
          std::string rev(end.rbegin(), end.rend());
          *current = rev;
        }
      } else {
        // Sequence does not contain a full k-mer. Need to recreate it
        // from the starting k-mer and seq by shifting.
        *current = start;
        for(auto it = seq.begin(); it != seq.end(); ++it)
          if(dir == forward) {
            current->shift_left(*it);
          } else {
            current->shift_right(*it);
          }
      }
    }

    // If the last k-mer has been used in a k-unitigs already, this
    // means two threads are working on the same unitigs, starting
    // from opposite ends. Output only if the current thread has the
    // "largest" end k-mer.
    bool end_new = insert_canonical(end_points_, *current);
    if(!end_new) {
      if(start.get_canonical() < current->get_canonical())
        return;
    }

    // Output results
    if(start.k() + seq.length() < args.min_len_arg)
      return;
    uint64_t id = (unitig_id_ += 1) - 1;
    output << ">" << id << " length:" << (start.k() + seq.size()) << "\n";
    if(dir == backward) {
      std::string reversed(seq.rbegin(), seq.rend());
      output << reversed << start << "\n";
    } else {
      output << start << seq << "\n";
    }
    output << jflib::endr;
  }
};
typedef bloom_filter<mer_type, mer_dna_hash, mt_access<unsigned int>> bloom_filter_type;
typedef create_k_unitig<mer_bloom_counter2, bloom_filter_type, bloom_filter_type, read_parser> unitiger_type;

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

  //  mer_type::k(args.mer_arg);
  if(!args.min_len_given)
    args.min_len_arg = args.mer_arg + 1;

  // Populate Bloom filter with k-mers
  std::auto_ptr<mer_bloom_counter2> kmers;

  {
    if(args.load_given) {
      mapped_file dbf(args.load_arg);
      uint64_t* base = (uint64_t*)dbf.base();
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
                           *output_ostream);
    unitiger.exec_join(args.threads_arg);
  }

  return EXIT_SUCCESS;
}
