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

#include <vector>

#include <jellyfish/dbg.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/mer_counting.hpp>

#include <src/error_correct_reads.hpp>

typedef uint64_t hkey_t;
typedef uint64_t hval_t;
//typedef jellyfish::invertible_hash::array<uint64_t,atomic::gcc,allocators::mmap> inv_hash_storage_t;
typedef std::vector<inv_hash_storage_t> hashes_t;

std::ostream &operator<<(std::ostream &os, const forward_counter &c) {
  return os << c._c;
}
std::ostream &operator<<(std::ostream &os, const backward_counter &c) {
  return os << c._c;
}



template<class instance_t>
class error_correct_t : public thread_exec {
  jellyfish::parse_read *_parser;
  hashes_t               *_hashes;
  int                     _mer_len;
  int                     _skip;
  int                     _good;
  int                     _anchor;
  std::string             _prefix;
  int                     _min_count;
  int                     _window;
  int                     _error;
public:
  error_correct_t(jellyfish::parse_read *parser, hashes_t *hashes) :
    _parser(parser), _hashes(hashes), 
    _mer_len(_hashes->front().get_key_len() / 2),
    _skip(0), _good(1), _min_count(1), _window(0), _error(0) {
  }

  void do_it(int nb_threads) {
    exec_join(nb_threads);
  }
  
  virtual void start(int id) {
    instance_t(this, id).start();
  }

  error_correct_t & skip(int s) { _skip = s; return *this; }
  error_correct_t & good(int g) { _good = g; return *this; }
  error_correct_t & anchor(int a) { _anchor = a; return *this; }
  error_correct_t & prefix(const char *s) { _prefix = s; return *this; }
  error_correct_t & prefix(const std::string s) { _prefix = s; return *this; }
  error_correct_t & mer_len(int l) { _mer_len = l; return *this; }
  error_correct_t & min_count(int m) { _min_count = m; return *this; }
  //  error_correct_t & advance(int a) { _advance = a; return *this; }
  error_correct_t & window(int w) { _window = w; return *this; }
  error_correct_t & error(int e) { _error = e; return *this; }

  inline hashes_t::const_iterator begin_hash() { return _hashes->begin(); }
  inline hashes_t::const_iterator end_hash() { return _hashes->end(); }

  jellyfish::parse_read *parser() const { return _parser; }
  int skip() const { return _skip; }
  int good() const { return _good; }
  int anchor() const { return _anchor; }
  const std::string & prefix() const { return _prefix; }
  int mer_len() const { return _mer_len; }
  int min_count() const { return _min_count; }
  //  int advance() const { return _advance; }
  int window() const { return _window ? _window : _mer_len; }
  int error() const { return _error ? _error : _mer_len / 2; }
};

class error_correct_instance {
public:
  typedef error_correct_t<error_correct_instance> ec_t ;

private:
  ec_t                     *_ec;
  int                       _id;
  size_t                    _buff_size;
  char                     *_buffer;
  hashes_t::const_iterator  chash;
  
public:
  error_correct_instance(ec_t *ec, int id) :
    _ec(ec), _id(id), _buff_size(0), _buffer(0) {}

  void start() {
    jellyfish::parse_read::thread parser = _ec->parser()->new_thread();
    const jellyfish::read_parser::read_t *read;
    std::ostringstream output_file;

    output_file << _ec->prefix() << "_" << _id;
    std::ofstream details((output_file.str() + ".log").c_str());
    if(!details.good())
      eraise(std::runtime_error)
        << "Failed to open file '" << output_file << ".log'" << err::no;
    std::ofstream output((output_file.str() + ".fa").c_str());
    if(!output.good())
      eraise(std::runtime_error)
        << "Failed opening file '" << output_file << ".fa'" << err::no;
    details.exceptions(std::ios::eofbit|std::ios::failbit|std::ios::badbit);
    output.exceptions(std::ios::eofbit|std::ios::failbit|std::ios::badbit);

    uint64_t nb_reads = 0;
    while((read = parser.next_read())) {
      nb_reads++;
      insure_length_buffer(read->seq_e - read->seq_s);
      
      kmer_t      mer;
      const char *input = read->seq_s + _ec->skip();
      char       *out   = _buffer + _ec->skip();
      DBG << V(_ec->skip()) << V((void*)read->seq_s) << V((void*)input);
      //Prime system. Find and write starting k-mer
      chash = _ec->begin_hash();
      if(!find_starting_mer(mer, input, read->seq_e, out)) {
        details << "Skipped " << substr(read->header, read->hlen) << "\n";
        continue;
      }
      DBG << V((void*)read->seq_s) << V((void*)input) << V(kmer_t::k());
      // Extend forward and backward
      forward_log fwd_log(_ec->window(), _ec->error());
      char *end_out = 
        extend(forward_mer(mer), forward_ptr<const char>(input),
               forward_counter(input - read->seq_s),
               forward_ptr<const char>(read->seq_e),
               forward_ptr<char>(out), fwd_log);
      DBG << V((void*)end_out) << V((void*)read->seq_e);
      assert(input > read->seq_s + kmer_t::k());
      assert(out > _buffer + kmer_t::k());
      assert(input - read->seq_s == out - _buffer);
      backward_log bwd_log(_ec->window(), _ec->error());
      char *start_out =
        extend(backward_mer(mer), 
               backward_ptr<const char>(input - kmer_t::k() - 1),
               backward_counter(input - kmer_t::k() - read->seq_s - 1),
               backward_ptr<const char>(read->seq_s - 1),
               backward_ptr<char>(out - kmer_t::k() - 1), bwd_log);
      DBG << V((void*)start_out) << V((void*)read->seq_s);
      start_out++;
      assert(start_out >= _buffer);
      assert(_buffer + _buff_size >= end_out);

      output << ">" << substr(read->header, read->hlen) 
             << " " << fwd_log << " " << bwd_log << "\n"
             << substr(start_out, end_out) << "\n";
    }
    details << "Nb reads " << nb_reads << std::endl;
    details.close();
    output.close();
  }

private:

  // Extend and correct read. Copy from input to out. mer should be
  // represent a "good" starting k-mer at the input position.
  // out point to the next character to be written.
  template<typename dir_mer, typename in_dir_ptr, typename out_dir_ptr,
           typename counter, typename elog>
  char * extend(dir_mer mer, in_dir_ptr input, 
                counter pos, in_dir_ptr end,
                out_dir_ptr out, elog &log) {
    counter cpos = pos;
    DBG << V((void*)input.ptr()) << V((void*)end.ptr()) << V(cpos);
    for( ; input < end; ++input) {
      char     base        = *input;
      DBG << V((void*)input.ptr()) << V((void*)end.ptr()) << V(base);
      if(base == '\n')
        continue;
      cpos = pos;
      ++pos;

      chash = _ec->begin_hash();
      uint64_t ori_code;
      if(!mer.shift(base)) {
        ori_code = 5; // Invalid base
        mer.shift((uint64_t)0);
      } else {
        ori_code = mer.code(0);
      }
      uint64_t counts[4];
      uint64_t ucode;
      int      count;
      
      while(true) {
        ucode = 0;
        count = get_all_alternatives(mer, counts, ucode);
        DBG << V(*cpos) << V(count) << V(counts[0]) << V(counts[1]) << V(counts[2]) << V(counts[3]);

        if(count == 0) {
          if(++chash == _ec->end_hash()) {
            log.truncation(cpos);
            goto done; // No continuation -> stop
          }
        } else
          break;
      }
      if(count == 1) { // One continuation. Is it an error?
        if(ucode != ori_code) {
          mer.replace(0, ucode);
          if(log.substitution(cpos, base, mer.base(0)))
            goto truncate;
        }
        *out++ = mer.base(0);
        continue;
      }

      // Check that there is at least one more base in the
      // sequence. If not, leave it along
      if(input >= end) {
        log.truncation(cpos);
        goto done;
      }

      // Select the replacement base to try. Find the one with highest
      // coverage and check that it is at least 3 times as big as
      // current base count and that it can continue one more base. In
      // that case, make the switch.
      uint64_t check_code = ori_code;
      uint64_t ori_count  = ori_code < 4 ? counts[ori_code] : 0;
      if(ori_count < (uint64_t)_ec->min_count())
        ori_count = 0;

      uint64_t max_count = 10000000000;
      if(ori_count <4)
        {
        max_count  = 3 * ori_count;
        }

      for(int i = 0; i < 4; i++) {
        if(counts[i] < (uint64_t)_ec->min_count())
            continue;
        if(counts[i] > max_count && counts[i]-ori_count<200) {
          check_code = i;
          max_count  = counts[i];
        }
      }
      if(check_code == ori_code) {
        // Don't need to check that check_code == 5 as an alternative
        // would have been found by now.
        *out++ = base;
        continue;
      }

      // Check that it continues at least one more base
      dir_mer    nmer   = mer;
      //      in_dir_ptr ninput = input;
      nmer.replace(0, check_code);
      // char       nbase  = *ninput++;
      // Does not matter what we shift, check all alternative anyway.
      nmer.shift((uint64_t)0);

      uint64_t   ncounts[4];
      int        ncount;
      uint64_t   nucode = 0;
      ncount = get_all_alternatives(nmer, ncounts, nucode);
      DBG << V(*cpos) << V(ncount) << V(ncounts[0]) << V(ncounts[1]) << V(ncounts[2]) << V(ncounts[3]);
      if(ncount > 0) {
        mer.replace(0, check_code);
        *out++ = mer.base(0);
        if(check_code != ori_code)
          if(log.substitution(cpos, base, mer.base(0)))
            goto truncate;
        // if(ncount == 1) { // While we are at it, there is a uniq continuation
        //   mer    = nmer;
        //   input  = ninput;
        //   ++pos;
        //   if(nucode != mer.code(0)) {
        //     mer.replace(0, nucode);
        //     if(log.substitution(cpos, nbase, mer.base(0)))
        //       goto truncate;
        //   }
        //   *out++ = mer.base(0);
        // }
      }
    }
    
  done:
    return out.ptr();

  truncate:
    int diff = log.remove_last_window();
    out = out - diff;
    DBG << V(*cpos) << V(diff) << V(*(cpos - diff));
    log.truncation(cpos - diff);
    goto done;
  }

  template <typename dir_mer>
  int get_all_alternatives(const dir_mer &mer, uint64_t counts[], uint64_t &ucode) {
    uint64_t val   = 0;
    dir_mer  nmer(mer);
    int      count = 0;
    for(uint64_t i = 0; i < (uint64_t)4; ++i) {
      nmer.replace(0, i);
      if((*chash).get_val(nmer.canonical(), val, true)) {
        counts[i] = val;
        if(val >= (uint64_t)_ec->min_count()) {
          count++;
          ucode = i;
        }
      } else {
        counts[i] = 0;
      }
    }
  
    return count;
  }

  void insure_length_buffer(size_t len) {
    if(len > _buff_size) {
      _buff_size = len > 2 * _buff_size ? len + 100 : 2 * _buff_size;
      _buffer = (char *)realloc(_buffer, _buff_size);

      if(!_buffer)
        eraise(std::runtime_error)
          << "Buffer allocation failed, size " << _buffer << err::no;
    }
  }

  bool find_starting_mer(kmer_t &mer, const char * &input, const char *end, char * &out) {
    while(input < end) {
      for(int i = 0; input < end && i < _ec->mer_len(); ++i) {
        char base = *input++;
        *out++ = base;
        if(!mer.shift_left(base))
          i = -1;        // If an N, skip to next k-mer
      }
      int found = 0;
      while(input < end) {
        hval_t val = 0;
        if(!(*chash).get_val(mer.canonical(), val, true))
          val = 0;
     
        found = (int)val >= _ec->anchor() ? found + 1 : 0;
        if(found >= _ec->good())
          return true;
        char base = *input++;
        *out++ = base;
        if(!mer.shift_left(base))
          break;
      }
    }
    return false;
  }
};

int main(int argc, char *argv[])
{
  args_t args(argc, argv);

  // Open Jellyfish databases
  hashes_t hashes;
  unsigned int key_len = 0;
  for(args_t::db_arg_const_it it = args.db_arg.begin(); it != args.db_arg.end(); ++it) {
    mapped_file dbf(*it);
    dbf.random().will_need();
    hashes.push_back(*raw_inv_hash_query_t(dbf).get_ary());
    
    if(key_len == 0)
      key_len = hashes.front().get_key_len();
    else if(key_len != hashes.back().get_key_len())
      die << "Different key length (" << hashes.back().get_key_len() 
          << " != " << key_len
          << ") for hash '" << *it << "'";
  }
  
  jellyfish::parse_read parser(args.file_arg.begin(), args.file_arg.end(), 100);

  kmer_t::k(key_len / 2);
  error_correct_instance::ec_t correct(&parser, &hashes);
  correct.skip(args.skip_arg).good(args.good_arg)
    .anchor(args.anchor_count_given ? args.anchor_count_arg : args.min_count_arg)
    .prefix(args.output_arg).min_count(args.min_count_arg)
    .window(args.window_given ? args.window_arg : key_len / 2)
    .error(args.error_given ? args.error_arg : key_len / 4);
  correct.do_it(args.thread_arg);

  return 0;
}
