#ifndef __READ_PARSER_HPP__
#define __READ_PARSER_HPP__

#include <fstream>
#include <memory>
#include <jflib/multiplexed_parser.hpp>
#include <charb.hpp>
//#include <err.hpp>

// Return reads from a fasta or fastq file. The parsing is done in a
// dedicated thread.

struct read_parser_read {
  charb header;
  charb sequence;
  charb quals;
};


class read_parser : public multiplexed_parser<read_parser_read> {
  typedef multiplexed_parser<read_parser_read> super;
  std::vector<std::filebuf*>      filebufs_;
  std::istream                    input_;
  bool                            close_input_;
  pthread_t                       reader_id_;
  const char*                     error_;

  std::filebuf* open_file(const char* path) {
    auto res = new std::filebuf();
    res->open(path, std::ios::in);
    std::string err("Failed to open file ");
    err += path;
    if(!res->is_open())
      throw std::runtime_error(err);
      //      eraise(std::runtime_error) << "Failed to open file '" << path << "'" << err::no;
    return res;
  }
public:
  typedef read_parser_read read;
  /** Read parser reading file path, with given expected number of
      threads (equivalently, number of conc_iterator used
      concurrently). Specifying too low a number of threads can
      results in poor performance (thread waiting on a lock
      constantly). group_size specify the number of reads grouped
      together. The default value is most likely sufficient.

      The total number of buffer created is 3 * nb_threads * group_size
   */
  explicit read_parser(const char* path, int nb_threads = 16, int group_size = 100) :
    super(nb_threads, group_size),
    input_(open_file(path)), close_input_(true)
  { start_parsing(); }

  /** Same as above reading from an already open istream. In this case
      the stream is not closed by this class destructor.
   */
  read_parser(std::istream& input, int nb_threads = 16, int group_size = 100) :
    super(nb_threads, group_size),
    input_(input.rdbuf()), close_input_(false)
  { start_parsing(); }

  /** Open multiple files
   */
  template<typename Iterator>
  read_parser(Iterator file_start, Iterator file_end, 
              int nb_threads = 16, int group_size = 100) :
    super(nb_threads, group_size), input_(0)
  {
    for( ; file_start != file_end; ++file_start)
      filebufs_.push_back(open_file(*file_start));
    start_parsing();
  }

  virtual ~read_parser();

private:
  // Start the appropriate reader loop based on examining the beginning of the file
  virtual void parser_loop();

  bool parse_input_stream();
  // Main loop parsing fasta & fastq format
  bool fasta_reader_loop();
  bool fastq_reader_loop();
};

#endif /* __READ_PARSER_HPP__ */
