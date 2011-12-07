/*  This file is part of k_unitig.

    k_unitig is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    k_unitig is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with k_unitig.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <jellyfish/err.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/allocators_mmap.hpp>
#include <jellyfish/thread_exec.hpp>
#include <signal.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>

#include <charbuf.hpp>
#include <src/create_k_unitigs_cmdline.hpp>
#include <jflib/pool.hpp>
#include <gzip_stream.hpp>
#include <src/aligned_simple_array.hpp>

// Structure and thread for writing out results
struct output_pair {
  charstream seq;
  charstream counts;
};
typedef jflib::pool<output_pair> output_pool;

typedef uint64_t hkey_t;
typedef uint64_t hval_t;

class key_count_array {
  struct kvp {
    hkey_t key;
    hval_t val;
  };
  size_t      index;
  size_t      capacity;
  struct kvp *storage;

  void enlarge() {
    capacity *= 2;
    storage = (struct kvp *)::realloc(storage, sizeof(struct kvp) * capacity);
  }
public:
  key_count_array() : index(0), capacity(16), storage(0) {
    storage = (struct kvp *)realloc(storage, sizeof(struct kvp) * capacity);
  }

  void push_back(const hkey_t k, const hval_t v) {
    const size_t i = index++;
    if(i >= capacity)
      enlarge();
    storage[i].key = k;
    storage[i].val = v;
  }

  void clear() { index = 0; }
  size_t size() const { return index; }
  hkey_t key(size_t i) { return storage[i].key; }
  hkey_t val(size_t i) { return storage[i].val; }
  void chop(size_t n) { index -= (n > index ? index : n); }
};


typedef std::pair<hkey_t,hval_t> mer_count;
//typedef std::vector<mer_count> seq_t;
typedef key_count_array seq_t;
//typedef jellyfish::compacted_hash::query<hkey_t,hval_t> hash_query_t;
// typedef jellyfish::invertible_hash::array<uint64_t,atomic::gcc,allocators::mmap> inv_hash_storage_t;
typedef inv_hash_storage_t::iterator iterator_t;
typedef simple_array::array<unsigned char,atomic::gcc,allocators::mmap> atomic_array;

class counter_t {
  volatile uint64_t count;
  atomic::gcc       atomic;

public:
  counter_t() : count(0) {}

  inline uint64_t operator++(int) {
    return atomic.fetch_add(&count, (uint64_t)1);
  }
  inline uint64_t inc(uint64_t x) {
    return atomic.fetch_add(&count, x);
  }
  inline uint64_t get() { return count; }

  class block {
    counter_t *c;
    uint64_t   bs;
    uint64_t   base, i;

    //    friend counter_t;
  public:
    block(counter_t *_c, uint64_t _bs) : c(_c), bs(_bs), base(0), i(bs) {}

  public:
    inline uint64_t operator++(int) {
      if(i >= bs) {
        i = 0;
        base = c->inc(bs);
      }
      return base + i++;
    }
  };
  block get_block(uint64_t bs = 100) { return block(this, bs); }
};

class collapse_kmers : public thread_exec {
  const create_k_unitigs_args  args;
  const inv_hash_storage_t    *hash;
  const hkey_t                 mask, forward_mask, backward_mask;
  const uint_t                 backward_shift;
  uint64_t                     uid_start;
  uint_t                       mer_len;
  atomic_array                 used_by;
  bool                         canonical;
  unsigned long                min_nb_kmers;
  counter_t                    slice_counter;
  counter_t                    fragment_counter;
  jflib::locks::barrier        done_barrier;

  std::ostream                *out_fasta;
  std::ostream                *out_counts;
  output_pool                  pool;

  const static int dir_backward = -1;
  const static int dir_forward  = 1;
  static collapse_kmers *this_progress;

  struct counters {
    uint64_t added, outputed, aborted, slices;
  };
  volatile struct counters *stats;
  int              nb_threads;
  size_t           nb_slices;

public:
  collapse_kmers(const inv_hash_storage_t *_hash, 
                 const create_k_unitigs_args &_args) :
    args(_args), hash(_hash),
    mask(((hkey_t)1 << hash->get_key_len()) - 1),
    forward_mask(mask - 3), backward_mask(mask >> 2),
    backward_shift(hash->get_key_len() - 2),
    mer_len(hash->get_key_len() / 2),
    used_by(hash->get_size()),
    canonical(args.both_strands_flag),
    min_nb_kmers(args.min_len_arg > mer_len ? args.min_len_arg + 1 - mer_len : 0),
    done_barrier(args.threads_arg),
    out_fasta(0), out_counts(0), pool(3 * args.threads_arg)
  {
    std::stringstream output_file_fasta, output_file_counts;
    output_file_fasta << args.prefix_arg << ".fa";
    output_file_counts << args.prefix_arg << ".counts";
    if(args.gzip_flag) {
      output_file_fasta << ".gz";
      output_file_counts << ".gz";
    }
    if(args.gzip_flag) {
      out_fasta = new gzipstream(output_file_fasta.str().c_str());
      if(args.counts_flag)
        out_counts = new gzipstream(output_file_counts.str().c_str());
    } else {
      out_fasta = new std::ofstream(output_file_fasta.str().c_str());
      if(args.counts_flag)
        out_counts = new std::ofstream(output_file_counts.str().c_str());
    }
    out_fasta->exceptions(std::ifstream::eofbit|std::ifstream::failbit|std::ifstream::badbit);
    if(args.counts_flag)
      out_counts->exceptions(std::ifstream::eofbit|std::ifstream::failbit|std::ifstream::badbit);
  }

  ~collapse_kmers() {
    delete out_fasta;
    delete out_counts;
  }

  static void alarm_handler(int sig) {
    if(this_progress) {
      this_progress->display_progress();
      alarm(2);
    }
  }

  void display_progress() {
    struct counters totals;
    memset(&totals, '\0', sizeof(totals));
    //volatile struct counters *nstats = stats;
    for(int i = 0; i < nb_threads; i++) {
      totals.added    += stats[i].added;
      totals.aborted  += stats[i].aborted;
      totals.outputed += stats[i].outputed;
      totals.slices   += stats[i].slices;
    }
    std::cerr << "\r slices " << totals.slices << " added " << totals.added << " output " << totals.outputed \
              << " abort " << totals.aborted << std::flush;
  }

  std::string mer_to_s(hkey_t key) const {
    char mer[33];
    jellyfish::parse_dna::mer_binary_to_string(key, mer_len, mer);
    return std::string(mer);
  }

  void do_it(int _nb_threads) {
    if(args.progress_flag) {
      this_progress = this;
      struct sigaction act;
      memset(&act, '\0', sizeof(act));
      act.sa_handler = alarm_handler;
      sigaction(SIGALRM, &act, NULL);
      alarm(2);
    }

    nb_threads = _nb_threads;
    nb_slices = nb_threads * 100;
    stats = new struct counters[nb_threads];
    if(args.progress_flag)
      std::cerr << "threads " << nb_threads << " slices " << nb_slices << std::endl;
    for(int i = 0; i < nb_threads; i++)
      memset((void *)&stats[i], '\0', sizeof(struct counters));
    exec_join(nb_threads + 1);

    if(args.progress_flag) {
      alarm(0);
      this_progress = 0;
      display_progress();
      std::cerr << "\n";
    }
    delete stats;
  }

  // What can happen if there is no unique continuation
  static const char * const dry; // no continuation
  static const char * const famb; // forward ambiguous
  static const char * const ramb;
  static const char * const fbranch; // forward branch
  static const char * const rbranch;
  static const char * const ful; // forward unique low count
  static const char * const rul;

  void output_thread() {
    while(true) {
      output_pool::elt e(pool.get_B());
      if(e.is_empty())
        break;
      out_fasta->write(e->seq.str(), e->seq.size());
      if(out_counts)
        out_counts->write(e->counts.str(), e->counts.size());
    }
  }


  // Start the threads. If th_id == nb_threads, this is the output thread.
  void start(int th_id) {
    if(th_id == nb_threads) {
      output_thread();
      return;
    }
      
    volatile struct counters *my_counters = &stats[th_id];

    counter_t::block frag_id = fragment_counter.get_block();
    const char *forward_term, *backward_term;
    seq_t forward, backward;
    uint64_t min_cov = args.min_cov_arg;

    for(size_t i = slice_counter++; i <= nb_slices; i = slice_counter++) {
      my_counters->slices++;
      iterator_t it = hash->iterator_slice(i, nb_slices);

      while(it.next()) {
        forward.clear();
        backward.clear();

        if(it.get_val() < min_cov) {
          my_counters->aborted++;
          continue;
        }
        forward.push_back(it.get_key(), it.get_val());
        forward_term = grow_sequence((unsigned char)(th_id & 0xff), it.get_key(),
                                     it.get_id(),
                                     &forward, dir_forward, &my_counters->added);
        if(!forward_term) {
          my_counters->aborted++;
          continue;
        }

        backward_term = grow_sequence((unsigned char)(th_id & 0xff), it.get_key(), 
                                      it.get_id(),
                                      &backward, dir_backward, &my_counters->added);
        if(!backward_term) {
          my_counters->aborted++;
          continue;
        }
        if((backward.size() + forward.size()) >= min_nb_kmers) {
          my_counters->outputed++;
          output_sequence(frag_id++, &forward, &backward,
                          forward_term, backward_term);
        }
      }
    }
    if(done_barrier.wait() == PTHREAD_BARRIER_SERIAL_THREAD)
      pool.close_A_to_B();
  }
  
  /* Return true if sequence has grown. Return (char *)0 if encounter a
     k-mer with id less than start_id (avoid exploring sequence
     already explored). Otherwise, return one of the error code explaining why it stopped.
   */
  const char * grow_sequence(unsigned char th_id, const hkey_t start_key, 
                             const uint64_t start_id,
                             seq_t *vec, const int direction, 
                             volatile uint64_t *added) {
    hkey_t      current_key = start_key;
    hkey_t      next_key = 0, cannon_next_key, bnext_key = 0, key_unused;
    hval_t      next_val = 0, val_unused;
    uint64_t    next_id = 0, id_unused;
    const char *found;
    size_t      low_cont = 0;
    size_t      low_stretch = args.low_stretch_arg;

    unsigned char current = used_by.set_if_empty(start_id, th_id + 1);
    if(direction == dir_forward) {
      // Get started only if first k-mer is empty
      if(current)
        return 0;
      (*added)++;
    } else {
      // Do backward direction only if starting k-mer still is 'ours'
      if(current != th_id + 1)
        return 0;
    }

    while(true) {
      found = unique_continuation(current_key, &next_key, &cannon_next_key, 
                                  &next_id, &next_val, direction);
      // TODO: This is a mess. Rethink this state checking and error reporting!
      if(found == 0) {
        low_cont = 0;
      } else if(found == dry) {
        if(low_cont)
          vec->chop(low_cont);
        return dry;
      } else if(found == ful) {
        if(args.cont_on_low_flag) {
          if(low_cont >= low_stretch) {
            vec->chop(low_cont);
            return dry;
          }
          ++low_cont;
        } else
          return dry;
      } else {
        return found;
      }
      
      const char *pfound = found;
      found = unique_continuation(next_key, &bnext_key, &key_unused,
                                  &id_unused, &val_unused, -direction);

      if(found == 0 || found == ful) {
        // good
      } else if(found == dry) {
        if(low_cont)
          vec->chop(low_cont);
        return dry;
      } else if(found == famb) {
        return ramb;
      } else if(found == fbranch) {
        return rbranch;
      } else {
        return found;
      }

      if(current_key != bnext_key) {
        if(args.cont_on_low_flag)
          return rbranch;
        char mer1_string[33], mer2_string[33], mer3_string[33];
        jellyfish::parse_dna::mer_binary_to_string(current_key, mer_len, mer1_string);
        jellyfish::parse_dna::mer_binary_to_string(bnext_key, mer_len, mer2_string);
        jellyfish::parse_dna::mer_binary_to_string(next_key, mer_len, mer3_string);
        std::cerr << "Backward unique next key not equal to current: " << mer1_string << " != " << mer2_string << " next was " << mer3_string << " pfound '" << (pfound ? pfound : "") << "' found '" << (found ? found : "") << "'" << " low_cont " << low_cont << std::endl;
        die << "Failed";
      }

      vec->push_back(next_key, next_val);
      current_key        = next_key;

      // If the next k-mer already belongs to another thread, trump it
      // if we have a larger thread id. Otherwise bail out.
      unsigned char omax = used_by.set_to_max(next_id, th_id + 1);
      if((th_id + 1) <= omax)
        return 0;
      if(!omax)
        (*added)++;
    }
  }

  // return (char *)0 if there is a unique continuation. Otherwise,
  // return one of the error code string.
  const char * unique_continuation(const hkey_t in_key, hkey_t *out_key, 
                                   hkey_t *cannon_out_key, uint64_t *out_id,
                                   hval_t *value, 
                                   const int direction)
  {
    int      found    = 0;
    int      lfound   = 0;
    int      lmfound  = 0;
    hval_t   res;
    uint64_t id;
    hkey_t   try_key, akey;
    uint64_t min_cov  = args.min_cov_arg;
    uint64_t min_cont = args.min_cont_arg;

    hkey_t   lc_key = 0, lc_cannon_key = 0;
    uint64_t lc_id = 0;
    hval_t   lc_val = 0;

    for(hkey_t i = 0; i < 4; i++) {
      if(direction == dir_forward)
        try_key = ((in_key << 2) & forward_mask) | i;
      else
        try_key = ((in_key >> 2) & backward_mask) | (i << backward_shift);

      if(canonical) {
        hkey_t rev = jellyfish::parse_dna::reverse_complement(try_key, mer_len);
        akey = (rev < try_key ? rev : try_key);
      } else {
        akey = try_key;
      }
      
      if(hash->get_val(akey, id, res, true)) {
        if(res >= min_cov) {
          ++lfound;
          if(res >= min_cont) {
            if(++found > 1)
              return fbranch;
            *out_key        = try_key;
            *cannon_out_key = akey;
            *out_id         = id;
            *value          = res;
          }
        } else {
          if(++lmfound == 1) {
            lc_key        = try_key;
            lc_cannon_key = akey;
            lc_id         = id;
            lc_val        = res;
          }
        }
      }
    }
    
    if(lfound > 1) return famb;
    if(found) return 0;
    if(lmfound == 1 && lfound == 0) {
      *out_key        = lc_key;
      *cannon_out_key = lc_cannon_key;
      *out_id         = lc_id;
      *value          = lc_val;
      return ful;
    }
    return dry;
  }

  void output_sequence(uint64_t id, seq_t *forward, seq_t *backward, 
                       const char *forward_term, const char *backward_term) {
    char mer_string[33];
    char trans[4] = {'A', 'C', 'G', 'T'};

    int fwd_start_id = 0;
    hkey_t first_mer;
    hval_t first_count;
    if(backward->size() > 0) {
      first_mer = backward->key(backward->size() - 1);
      first_count = backward->val(backward->size() - 1);
    } else {
      first_mer = forward->key(0);
      first_count = forward->val(0);
      ++fwd_start_id;
    }

    // Compute the sum of the qual values
    uint64_t sumq = first_count;
    for(int i = backward->size() - 2; i >= 0; --i)
      sumq += backward->val(i);
    for(uint_t i = fwd_start_id; i < forward->size(); ++i)
      sumq += forward->val(i);

    mer_string[32] = '\0';
    
    output_pool::elt e(pool.get_A());
    e->seq.rewind();
    e->counts.rewind();
    (e->seq) << ">" <<  id 
           << " length:" << (mer_len - 1 + backward->size() + forward->size()) 
           << " bwd:" << backward_term << " fwd:" << forward_term 
           << " sumq:" << sumq << "\n";
    e->counts << ">" << id << "\n";

    jellyfish::parse_dna::mer_binary_to_string(first_mer, mer_len, mer_string);
    e->seq << mer_string;

    for(uint_t i = 1; i < mer_len; ++i)
      e->counts << ". ";
    e->counts << first_count << " ";
    size_t total = mer_len + 1;
    for(int i = backward->size() - 2; i >= 0; --i, ++total) {
      e->seq << trans[backward->key(i) & 0x3];
      e->counts << backward->val(i) << " ";
      if(total % 70 == 0) {
        e->seq << "\n";
        e->counts << "\n";
      }
    }
    for(uint_t i = fwd_start_id; i < forward->size(); ++i, ++total) {
      e->seq << trans[forward->key(i) & 0x3];
      e->counts << forward->val(i) << " ";
      if(total % 70 == 0) {
        e->seq << "\n";
        e->counts << "\n";
      }
    }
    if(total % 70 != 1) {
      e->seq << "\n";
      e->counts << "\n";
    }
  }
};
collapse_kmers *collapse_kmers::this_progress = 0;
const char * const collapse_kmers::dry = "dry";
const char * const collapse_kmers::famb = "ambiguous_branch";
const char * const collapse_kmers::ramb = "ambiguous_merge";
const char * const collapse_kmers::fbranch = "branch";
const char * const collapse_kmers::rbranch = "merge";
const char * const collapse_kmers::ful = "unique_low";
const char * const collapse_kmers::rul = "unique_low";


int main(int argc, char *argv[])
{
  create_k_unitigs_args args(argc, argv);

  mapped_file dbf(args.file_arg);
  dbf.random().will_need().load();
  raw_inv_hash_query_t hash(dbf);

  // By default, the minimum length is k+1
  if(!args.min_len_given)
    args.min_len_arg = hash.get_mer_len() + 1;

  collapse_kmers(hash.get_ary(), args).do_it(args.threads_arg);
  return 0;
}
