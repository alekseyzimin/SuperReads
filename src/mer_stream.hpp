#ifndef __MER_STREAM_HPP__
#define __MER_STREAM_HPP__

template<typename mer_type, typename parser_type>
class mer_stream {
  mer_type                     fmer_;
  mer_type                     rmer_;
  typename parser_type::stream input_stream_;
  unsigned int                 bases_;
  const char*                  s_;
  const char*                  e_;

public:
  mer_stream(int mer_len, parser_type& parser) :
    fmer_(mer_len), rmer_(mer_len), input_stream_(parser),
    bases_(0), s_(0), e_(0)
  {
    if(input_stream_) {
      s_ = input_stream_->sequence;
      e_ = input_stream_->sequence.end();
      this->operator++();
    }
  }

  mer_stream& operator++() {
    while(true) {
      while(s_ < e_) {
        uint64_t code = mer_type::code(*s_++);
        switch(code) {
        case mer_type::CODE_COMMENT:
          s_ = e_; // Skip read
          break;
        case mer_type::CODE_RESET:
          bases_ = 0;
          break;
        case mer_type::CODE_IGNORE:
          break;
        default:
          fmer_.shift_left(code);
          rmer_.shift_right(mer_type::complement(code));
          if(++bases_ >= fmer_.k())
            return *this;
          break;
        }
      }

      ++input_stream_;
      if(!input_stream_) {
        s_ = e_ = 0;
        bases_ = 0;
        return *this;
      }
      s_     = input_stream_->sequence;
      e_     = input_stream_->sequence.end();
      bases_ = 0;
    }
  }

  operator void*() const { return (void*)s_; }

  const mer_type& operator*() const { return fmer_; }
  const mer_type* operator->() const { return &fmer_; }
  const mer_type& fmer() const { return fmer_; }
  const mer_type& rmer() const { return rmer_; }
  const mer_type& canonical() const { return fmer_ < rmer_ ? fmer_ : rmer_; }
};


#endif /* __MER_STREAM_HPP__ */
