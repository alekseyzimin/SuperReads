#include <src/read_parser.hpp>

read_parser::~read_parser() {
  if(reader_started_) {
    pool_.close_B_to_A();
    pthread_cancel(reader_id_);
    pthread_join(reader_id_, 0);
  }
  if(close_input_)
    delete input_.rdbuf();
}

void read_parser::reader_loop() {
  switch(input_.peek()) {
  case '>': fasta_reader_loop(); break;
    //  case '@': fastq_reader_loop(); break;
  default:
    break; // Should never be reached
  }
}

void read_parser::start_parsing_thread() {
  // Check input format
  switch(input_.peek()) {
  case EOF: pool_.close_A_to_B(); return; // Empty file -> don't do anything
  case '>': break;
  case '@': eraise(std::runtime_error) << "Fastq parsing not yet implemented";
  default:
    eraise(std::runtime_error) << "Invalid input file format";
  }

  // Finish initialization of the read_groups
  for(auto it = pool_.begin(); it != pool_.end(); ++it) {
    it->reads.resize(group_size_);
    it->nb_filled = 0;
  }

  auto self  = new self_pointer;
  self->self = this;
  int  res   = pthread_create(&reader_id_, 0, start_reader_loop , (void*)self);
  if(res) {
    delete self;
    eraise(std::runtime_error) << "Failed to create the reader thread: "
                               << err::str(res);
  }
  reader_started_ = true;
}

void read_parser::fasta_reader_loop() {
  int nextc = input_.peek();
  while(nextc != EOF) {
    read_pool::elt e(pool_.get_A());
    if(e.is_empty())
      break;

    int i = 0;
    for( ; nextc != EOF && i < group_size_; ++i) {
      read& r = e->reads[i];
      getline(input_, r.header);
      r.sequence.clear();
      nextc = input_.peek();
      while(nextc != EOF && nextc != '>') {
        getline_append(input_, r.sequence);
        r.sequence.chomp();
        nextc = input_.peek();
      }
    }
    e->nb_filled = i;
  }
  pool_.close_A_to_B();

    // Code with read groups
    // int& i = e.nb_filled;
    // for(i = 0; i < group_size_; ++i) {
    //   read& r = e.reads[i];
    //   getline(input_, r.header);
    //   r.sequence.clear();
    //   int nextc = input_.peek();
    //   while(nextc != EOF && nextc != '>') {
    //     getline_append(input_, r.sequence);
    //     r.sequence.chomp();
    //     nextc = input_.peek();
    //   }
    // }
}

