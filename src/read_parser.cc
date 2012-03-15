#include <src/read_parser.hpp>

read_parser::~read_parser() {
  if(close_input_)
    delete input_.rdbuf();
}

void read_parser::parser_loop() {
  switch(input_.peek()) {
  case '>': fasta_reader_loop(); break;
  case '@': fastq_reader_loop(); break;
  default:
    super::report_error("Invalid file format");
    break; // Should never be reached
  }
}

void read_parser::fasta_reader_loop() {
  int nextc = input_.peek();
  while(nextc != EOF) {
    elt e(elt_init());
    if(e.is_empty())
      break;

    size_type& i = e->nb_filled;
    for(i = 0 ; nextc != EOF && i < group_size(); ++i) {
      read& r = e->elements[i];
      getline(input_, r.header);
      r.sequence.clear();
      nextc = input_.peek();
      while(nextc != EOF && nextc != '>') {
        getline_append(input_, r.sequence);
        r.sequence.chomp();
        nextc = input_.peek();
      }
    }
  }
}

void read_parser::fastq_reader_loop() {
  charb unused_line;
  int nextc = input_.peek();
  while(nextc != EOF) {
    elt e(elt_init());
    if(e.is_empty())
      break;

    size_type& i = e->nb_filled;
    for(i= 0; nextc != EOF && i < group_size(); ++i) {
      read& r = e->elements[i];
      getline(input_, r.header);
      if(r.header[0] != '@') {
        report_error("Found bad sequence header");
        nextc = EOF;
        break;
      }
      r.sequence.clear();
      r.quals.clear();
      nextc = input_.peek();
      while(nextc != EOF && nextc != '+') {
        getline_append(input_, r.sequence);
        r.sequence.chomp();
        nextc = input_.peek();
      }
      getline(input_, unused_line); // read quals header: ignored
      while(nextc != EOF && r.quals.size() < r.sequence.size()) {
        getline_append(input_, r.quals);
        r.quals.chomp();
        nextc = input_.peek();
      }
      if(r.quals.size() != r.sequence.size()) { // Invalid input file
        report_error("Number of qual values != number of bases");
        nextc = EOF;
        break;
      }
    }
  }
}

