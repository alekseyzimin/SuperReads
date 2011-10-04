#ifndef __ARGS_T_HPP__
#define __ARGS_T_HPP__

#include <yaggo.hpp>

class args_t {
public:
  std::vector<const char *>      db_arg;
  typedef std::vector<const char *>::iterator db_arg_it;
  typedef std::vector<const char *>::const_iterator db_arg_const_it;
  bool                           db_given;
  uint32_t                       thread_arg;
  bool                           thread_given;
  bool                           both_strands_flag;
  uint32_t                       min_count_arg;
  bool                           min_count_given;
  uint32_t                       skip_arg;
  bool                           skip_given;
  uint32_t                       good_arg;
  bool                           good_given;
  uint32_t                       anchor_count_arg;
  bool                           anchor_count_given;
  uint32_t                       window_arg;
  bool                           window_given;
  uint32_t                       error_arg;
  bool                           error_given;
  yaggo::string                  output_arg;
  bool                           output_given;
  std::vector<const char *>      file_arg;
  typedef std::vector<const char *>::iterator file_arg_it;
  typedef std::vector<const char *>::const_iterator file_arg_const_it;

  enum {
    USAGE_OPT = 1000
  };

  args_t(int argc, char *argv[]) :
    db_arg(), db_given(false),
    thread_arg(1), thread_given(false),
    both_strands_flag(true),
    min_count_arg(2), min_count_given(false),
    skip_arg(2), skip_given(false),
    good_arg(2), good_given(false),
    anchor_count_arg(0), anchor_count_given(false),
    window_arg(0), window_given(false),
    error_arg(5), error_given(false),
    output_arg("error_corrected"), output_given(false)
  {
    static struct option long_options[] = {
      {"db", 1, 0, 'd'},
      {"thread", 1, 0, 't'},
      {"both-strands", 0, 0, 'C'},
      {"min-count", 1, 0, 'm'},
      {"skip", 1, 0, 's'},
      {"good", 1, 0, 'g'},
      {"anchor-count", 1, 0, 'a'},
      {"window", 1, 0, 'w'},
      {"error", 1, 0, 'e'},
      {"output", 1, 0, 'o'},
      {"help", 0, 0, 'h'},
      {"usage", 0, 0, USAGE_OPT},
      {"version", 0, 0, 'V'},
      {0, 0, 0, 0}
    };
    static const char *short_options = "hVd:t:Cm:s:g:a:w:e:o:";

    std::string err;
#define CHECK_ERR(type,val,which) if(!err.empty()) { std::cerr << "Invalid " #type " '" << val << "' for [" which "]: " << err << "\n"; exit(1); }
    while(true) { 
      int index = -1;
      int c = getopt_long(argc, argv, short_options, long_options, &index);
      if(c == -1) break;
      switch(c) {
      case ':': 
        std::cerr << "Missing required argument for "
                  << (index == -1 ? std::string(1, (char)optopt) : std::string(long_options[index].name))
                  << std::endl;
        exit(1);
      case 'h':
        std::cout << usage() << "\n\n" << help() << std::endl;
        exit(0);
      case USAGE_OPT:
        std::cout << usage() << "\nUse --help for more information." << std::endl;
        exit(0);
      case 'V':
        print_version();
        exit(0);
      case '?':
        std::cerr << "Use --usage or --help for some help\n";
        exit(1);
      case 'd':
        db_given = true;
        db_arg.push_back(optarg);
        break;
      case 't':
        thread_given = true;
        thread_arg = yaggo::conv_uint<uint32_t>((const char *)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-t, --thread=uint32")
        break;
      case 'C':
        both_strands_flag = true;
        break;
      case 'm':
        min_count_given = true;
        min_count_arg = yaggo::conv_uint<uint32_t>((const char *)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-m, --min-count=uint32")
        break;
      case 's':
        skip_given = true;
        skip_arg = yaggo::conv_uint<uint32_t>((const char *)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-s, --skip=uint32")
        break;
      case 'g':
        good_given = true;
        good_arg = yaggo::conv_uint<uint32_t>((const char *)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-g, --good=uint32")
        break;
      case 'a':
        anchor_count_given = true;
        anchor_count_arg = yaggo::conv_uint<uint32_t>((const char *)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-a, --anchor-count=uint32")
        break;
      case 'w':
        window_given = true;
        window_arg = yaggo::conv_uint<uint32_t>((const char *)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-w, --window=uint32")
        break;
      case 'e':
        error_given = true;
        error_arg = yaggo::conv_uint<uint32_t>((const char *)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-e, --error=uint32")
        break;
      case 'o':
        output_given = true;
        output_arg.assign(optarg);
        break;
      }
    }
    if(!db_given)
      error("[-d, --db=jellyfish.db] required switch");
    if(argc - optind < 1)
      error("Requires at least 1 argument.");
    for( ; optind < argc; ++optind) {
      file_arg.push_back(argv[optind]);
    }
  }
#define args_t_USAGE "Usage: error_correct_reads [options] file:path+"
  const char * usage() const { return args_t_USAGE; }
  void error(const char *msg) { 
    std::cerr << "Error: " << msg << "\n" << usage()
              << "\nUse --help for more information"
              << std::endl;
    exit(1);
  }
#define args_t_HELP "Error correct reads from a fastq file based on the k-mer frequencies.\n\n" \
  "Options (default value in (), *required):\n" \
  " -d, --db=jellyfish.db                   *Jellyfish database\n" \
  " -t, --thread=uint32                      Number of threads (1)\n" \
  " -C, --both-strands                       Canonical k-mers in database (true)\n" \
  " -m, --min-count=uint32                   Minimum count for a k-mer to be considered \"good\" (2)\n" \
  " -s, --skip=uint32                        Number of bases to skip for start k-mer (2)\n" \
  " -g, --good=uint32                        Number of good k-mer in a row for anchor (2)\n" \
  " -a, --anchor-count=uint32                Minimum count for an anchor k-mer (default=min-count)\n" \
  " -w, --window=uint32                      Size of window (default=mer length)\n" \
  " -e, --error=uint32                       Maximum number of error in a window (5)\n" \
  " -o, --output=prefix                      Output file prefix (error_corrected)\n" \
  "     --usage                              Usage\n" \
  " -h, --help                               This message\n" \
  " -V, --version                            Version"

  const char * help() const { return args_t_HELP; }
#define args_t_HIDDEN "Hidden options:"

  const char * hidden() const { return args_t_HIDDEN; }
  void print_version(std::ostream &os = std::cout) const {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
    os << PACKAGE_VERSION << "\n";
  }
  void dump(std::ostream &os = std::cout) {
    os << "db_given:" << db_given << " db_arg:" << yaggo::vec_str(db_arg) << "\n";
    os << "thread_given:" << thread_given << " thread_arg:" << thread_arg << "\n";
    os << "both_strands_flag:" << both_strands_flag << "\n";
    os << "min_count_given:" << min_count_given << " min_count_arg:" << min_count_arg << "\n";
    os << "skip_given:" << skip_given << " skip_arg:" << skip_arg << "\n";
    os << "good_given:" << good_given << " good_arg:" << good_arg << "\n";
    os << "anchor_count_given:" << anchor_count_given << " anchor_count_arg:" << anchor_count_arg << "\n";
    os << "window_given:" << window_given << " window_arg:" << window_arg << "\n";
    os << "error_given:" << error_given << " error_arg:" << error_arg << "\n";
    os << "output_given:" << output_given << " output_arg:" << output_arg << "\n";
    os << "file_arg:" << yaggo::vec_str(file_arg) << "\n";
  }
private:
};

#endif // __ARGS_T_HPP__"
