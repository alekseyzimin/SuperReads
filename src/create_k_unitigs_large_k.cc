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

class mer_stream {
  mer_dna             fmer_, rmer_;
  read_parser::stream read_stream_;
  unsigned int        bases_;
  const char*         s_;
  const char*         e_;

public:
  mer_stream(int mer_len, read_parser& parser) :
    fmer_(mer_len), rmer_(mer_len), read_stream_(parser),
    bases_(0), s_(0), e_(0)
  {
    if(read_stream_) {
      s_ = read_stream_->sequence;
      e_ = read_stream_->sequence.end();
    }
  }

  mer_stream& operator++() {
    while(true) {
      for( ; s_ != e_; ++s_) {
        uint64_t code = mer_dna::code(*s_);
        if(code != mer_dna::bad_code) {
          fmer_.shift_left(code);
          rmer_.shift_right(mer_dna::complement(code));
          if(++bases_ >= fmer_.k())
            return *this;
        } else
          bases_ = 0;
      }

      ++read_stream_;
      if(!read_stream_) {
        s_ = e_ = 0;
        return *this;
      }
      s_ = read_stream_->sequence;
      e_ = read_stream_->sequence.end();
    }
  }

  operator void*() const { return (void*)s_; }

  const mer_dna& operator*() const { return fmer_; }
  const mer_dna* operator->() const { return &fmer_; }
  const mer_dna& fmer() const { return fmer_; }
  const mer_dna& rmer() const { return rmer_; }
  const mer_dna& canonical() const { return fmer_ < rmer_ ? fmer_ : rmer_; }
};

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
    mer_stream k_mers(mer_len_, parser_);

    for( ; k_mers; ++k_mers)
      ++set_[k_mers.canonical()];
  }
};
typedef populate_mer_set<mer_bloom_counter2> mer_populate;

/* - map_type is maps k-mer to counts (like a counting bloom
 *   counter). Has method operator[] returning an unsigned int.
 *
 * - used_type and printed_type are sets of k-mer (e.g. bloom filter
 *   or hash table).
 *
 * Every type is assumed to be thread safe
 */
template<typename map_type, typename used_type, typename printed_type>
class create_k_unitigs {
  int          mer_len_;
  map_type     counts_;
  used_type    used_;
  printed_type printed_;
  read_parser  parser_;

  static const int forward  = 1;
  static const int backward = -1;

public:
  create_k_unitigs(int mer_len, const map_type &counts, used_type& used, printed_type& printed,
                   int threads, const char* sequence_file, std::ostream& output) :
    counts_(counts), used_(used), printed_(printed), parser_(sequence_file, threads) { }
  virtual ~create_k_unitigs() { }

  void start(int thid) {
    mer_stream k_mers(mer_len_, parser_);
    mer_dna conts[4];

    for( ; k_mers; ++k_mers) {
      if(used_.insert(k_mers.canonical()))
        continue; // Already seen -> skip

      int fwd_continuations = unique_continuation(*k_mers, forward, conts);
      for(int i = 0; i < fwd_continuations; ++i)
        grow_k_unitig(conts[i], forward);

      int bwd_continuations = unique_continuation(*k_mers, backward, conts);
      for(int i = 0; i < bwd_continuations; ++i)
        grow_k_unitig(conts[i], backward);

      if(fwd_continuations == 0 ||
         (fwd_continuations > 1 && bwd_continuations == 1))
        grow_k_unitig(*k_mers, backward);

      if(bwd_continuations == 0 ||
         (bwd_continuations > 1 && fwd_continuations == 1))
        grow_k_unitig(*k_mers, forward);
    }
  }

private:
  // Find number of continuation of <m> in <direction> and return in
  // the array conts the possible continuations.
  int find_continuations(const mer_dna& m, int direction, mer_dna conts[]) {
    int     pos_base     = direction == forward ? 0 : m.k() - 1;
    int     rev_pos_base = direction == forward ? m.k() - 1 : 0;
    int     i            = 0;
    mer_dna revc;

    conts[0] = m
    if(direction == forward)
      conts[0].shift_left((uint64_t)0);
    else
      conts[0].shift_right((uint64_t)0);

    static const char bases[4] = { 'A', 'C', 'G', 'T' };
    for(auto it = bases; it != bases + 4; ++it) {
      conts[i].base(pos_base)     = *it;
      revc = conts[i];
      revc.reverse_complement();
      if(counts_[conts[i] < revc ? conts[i] : revc] > 1) {
        conts[i+1] = conts[i];
        ++i;
      }
    }
    return i;
  }

  void grow_k_unitig(const mer_dna& start_m, int direction) {
    if(printed_.insert(start_m.canonical()))
      return;

    std::vector<char> seq;
    mer_dna           m(start_m);
    mer_dna           conts[4];
    const char*       reason   = 0;
    int               pos_base = direction == forward ? 0 : m.k() - 1;


    while(true) {
      int continuations = unique_continuation(m, direction, conts);
      switch(continuations) {
      case 0: reason = dry; goto done;
      case 1: break;
      default: reason = fbranch; goto done;
      }

      m = conts[0];
      continuations = unique_continuation(m, -direction, conts);
        
      switch(continuations) {
      case 1: break;
      default: reason = fbranch; goto done;
      }

      seq.push_back(m.base(pos_base));
    }

  done:
    if(seq.empty())
      return;
    if(printed_.insert(m.canonical()))
    print_k_unitig(m, seq, reason);
  }

  void print_k_unitig(const mer_dna& m, std::vector<char>& seq, const char* reason) {
    
  }
};

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
