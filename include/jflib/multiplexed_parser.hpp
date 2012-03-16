#ifndef __MULTIPLEXED_PARSER_HPP__
#define __MULTIPLEXED_PARSER_HPP__

#include <pthread.h>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <jflib/pool.hpp>
// #include <err.hpp>

template <typename T>
class multiplexed_parser {
  struct group {
    typedef std::vector<T>                 elt_vector;
    typedef typename elt_vector::size_type size_type;
    
    elt_vector elements;
    size_type  nb_filled;
  };

  typedef jflib::pool<group>        pool;
public:
  typedef typename group::size_type size_type;

private:
  const size_type group_size_;
  pool            pool_;
  bool            parser_started_;
  pthread_t       parser_id_;
  const char*     error_;

public:
  explicit multiplexed_parser(int nb_threads = 16,  size_type group_size = 100) :
    group_size_(group_size), pool_(3 * nb_threads), 
    parser_started_(false), error_(0) 
  { }
  virtual ~multiplexed_parser();

  // good if no error and not end_of_file
  bool good() const { return !pool_.is_closed_A_to_B() && error() == 0; }
  bool eof() const { return pool_.is_closed_A_to_B(); }
  // fail if an error occurred. end_of_file does not set fail
  bool fail() const { return error() != 0; }
  // error message
  const char* error() const { return jflib::a_load_ptr((const char*&)error_); }

  // Stream of element
  typedef typename pool::elt elt;
  class stream { 
    pool&     pool_;
    elt elt_;
    size_type i_;
  public:
    explicit stream(multiplexed_parser& rp) :
      pool_(rp.pool_), elt_(pool_.get_B()), i_(0) 
    { 
      while(!elt_.is_empty() && elt_->nb_filled == (size_type)0)
        elt_ = pool_.get_B();
    }
    // Probably useless
    stream() = default;
    // Non copyable
    stream(const stream& rhs) = delete;
    stream& operator=(const stream& rhs) = delete;

    T& operator*() { return elt_->elements[i_]; }
    T* operator->() { return &elt_->elements[i_]; }

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

  void start_parsing();

  size_type group_size() const { return group_size_; }
  
protected:
  // Actual parsing
  virtual void parser_loop() = 0;

  void report_error(const char* msg) { jflib::a_store_ptr(error_, msg); }
  typedef typename pool::side side;
  side& elt_init() { return pool_.get_A(); }

private:
  struct self_pointer {
    multiplexed_parser* self;
  };
  static void* static_start_parser_loop(void* self_);
  void start_parser_loop();
};

template<typename T>
multiplexed_parser<T>::~multiplexed_parser() {
  if(parser_started_) {
    pool_.close_B_to_A();
    // XXX: Do we need to cancel the thread? This seems to leak memory and
    // may not be necessary or desirable.
    //    pthread_cancel(parser_id_);
    pthread_join(parser_id_, 0);
  }
}

template<typename T>
void* multiplexed_parser<T>::static_start_parser_loop(void* self_) {
  std::auto_ptr<self_pointer> self((self_pointer*)self_);
  self->self->start_parser_loop();
  return 0;
}

template<typename T>
void multiplexed_parser<T>::start_parser_loop() {
  parser_loop(); // Call virtual subclass method
  pool_.close_A_to_B();
}

template<typename T>
void multiplexed_parser<T>::start_parsing() {
  // Finish initialization of the read_groups
  for(auto it = pool_.begin(); it != pool_.end(); ++it) {
    it->elements.resize(group_size_);
    it->nb_filled = 0;
  }

  auto self  = new self_pointer;
  self->self = this;
  int  res   = pthread_create(&parser_id_, 0, static_start_parser_loop , (void*)self);
  if(res) {
    delete self;
    throw std::runtime_error("Failed to create the parser thread");
  }
  parser_started_ = true;
}

#endif /* __MULTIPLEXED_PARSER_HPP__ */
