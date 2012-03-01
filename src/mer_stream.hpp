#ifndef __MER_STREAM_HPP__
#define __MER_STREAM_HPP__

template<typename mer_type>
class mer_stream {
  mer_type            mer_;
  read_parser::stream read_stream_;
  unsigned int        bases_;
  const char*         s_;
  const char*         e_;

public:
  mer_stream(int mer_len, read_parser& parser) :
    mer_(mer_len), read_stream_(parser),
    bases_(0), s_(0), e_(0)
  {
    if(read_stream_) {
      s_ = read_stream_->sequence;
      e_ = read_stream_->sequence.end();
    }
  }
  // Unfinished
  mer_stream& operator++() {
    while(true) {
      for( ; s_ != e_; ++s_) {
        uint64_t code = mer_type::code(*s_);
        if(code != mer_type::bad_code) {
          mer_.shift_left(code);
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
};


#endif /* __MER_STREAM_HPP__ */
