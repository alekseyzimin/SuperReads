#ifndef __READ_PARSER_HPP__
#define __READ_PARSER_HPP__

#include <pthread.h>
#include <fstream>
#include <memory>
#include <jflib/pool.hpp>
#include <charb.hpp>
#include <err.hpp>

// Return reads from a fasta or fastq file. The parsing is done in a
// dedicated thread.

class read_parser {
  struct read {
    charb header;
    charb sequence;
  };
  struct read_group {
    std::vector<read> reads;
    int nb_filled;
  };
  typedef jflib::pool<read_group> read_pool;
  std::istream                    input_;
  bool                            close_input_;
  int                             group_size_;
  read_pool                       pool_;
  bool                            reader_started_;
  pthread_t                       reader_id_;

  std::filebuf* open_file(const char* path) {
    auto res = new std::filebuf();
    res->open(path, std::ios::in);
    if(!res->is_open())
      eraise(std::runtime_error) << "Failed to open file '" << path << "'" << err::no;
    return res;
  }
public:
  /** Read parser reading file path, with given expected number of
      threads (equivalently, number of conc_iterator used
      concurrently). Specifying too low a number of threads can
      results in poor performance (thread waiting on a lock
      constantly). group_size specify the number of reads grouped
      together. The default value is most likely sufficient.

      The total number of buffer created is 3 * nb_threads * group_size
   */
  read_parser(const char* path, int nb_threads = 16, int group_size = 100) :
    input_(open_file(path)), close_input_(true), group_size_(group_size),
    pool_(3 * nb_threads), reader_started_(false)
  { start_parsing_thread(); }

  /** Same as above reading from an already open istream. In this case
      the stream is not closed by this class destructor.
   */
  read_parser(std::istream& input, int nb_threads = 16, int group_size = 100) :
    input_(input.rdbuf()), close_input_(false), group_size_(group_size),
    pool_(3 * nb_threads), reader_started_(false)
  { start_parsing_thread(); }

  virtual ~read_parser();

  // Stream of reads
  class stream { 
    read_pool&     pool_;
    read_pool::elt elt_;
    int            i_;
  public:
    stream(read_parser& rp) : pool_(rp.pool_), elt_(pool_.get_B()), i_(0) { 
      while(!elt_.is_empty() && elt_->nb_filled == 0)
        elt_ = pool_.get_B();
    }
    // Probably useless
    stream() = default;
    // Non copyable
    stream(const stream& rhs) = delete;
    stream& operator=(const stream& rhs) = delete;

    read& operator*() { return elt_->reads[i_]; }
    read* operator->() { return &elt_->reads[i_]; }

    stream& operator++() {
      if(++i_ < elt_->nb_filled)
        return *this;
      i_ = 0;
      do {
        elt_ = pool_.get_B();
      } while(!elt_.is_empty() && elt_->nb_filled == 0);
      return *this;
    }
    operator void*() const { return elt_.is_empty() ? (void*)0 : (void*)&elt_; }
  };

private:
  struct self_pointer {
    read_parser* self;
  };
  static void* start_reader_loop(void* self_) {
    std::auto_ptr<self_pointer> self((self_pointer*)self_);
    self->self->reader_loop();
    return 0;
  }

  // Decide which format the file is and throw an error if not
  // support. Finish initialization and start the parsing thread.
  void start_parsing_thread();

  // Start the approriate reader loop based on examining the beginning of the file
  void reader_loop();

  // Main loop parsing fasta format
  void fasta_reader_loop();
};

#endif /* __READ_PARSER_HPP__ */
