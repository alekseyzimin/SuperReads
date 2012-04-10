/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/***** This code was generated by Yaggo. Do not edit ******/

#ifndef __JOINKUNITIGS_V3_HPP__
#define __JOINKUNITIGS_V3_HPP__

#include <yaggo.hpp>

class joinKUnitigs_v3 {
public:
  int                            min_overlap_length_arg;
  bool                           min_overlap_length_given;
  const char *                   mean_and_stdev_by_prefix_file_arg;
  bool                           mean_and_stdev_by_prefix_file_given;
  const char *                   unitig_lengths_file_arg;
  bool                           unitig_lengths_file_given;
  const char *                   kunitigs_translation_file_arg;
  bool                           kunitigs_translation_file_given;
  const char *                   overlaps_file_arg;
  bool                           overlaps_file_given;
  const char *                   num_kunitigs_file_arg;
  bool                           num_kunitigs_file_given;
  int                            max_nodes_allowed_arg;
  bool                           max_nodes_allowed_given;
  const char *                   output_arg;
  bool                           output_given;
  int                            threads_arg;
  bool                           threads_given;
  const char *                   input_file_arg;

  enum {
    USAGE_OPT = 1000,
    MIN_OVERLAP_LENGTH_OPT,
    KUNITIGS_TRANSLATION_FILE_OPT,
    NUM_KUNITIGS_FILE_OPT,
    MAX_NODES_ALLOWED_OPT
  };

  joinKUnitigs_v3() : 
    min_overlap_length_arg(), min_overlap_length_given(false),
    mean_and_stdev_by_prefix_file_arg(""), mean_and_stdev_by_prefix_file_given(false),
    unitig_lengths_file_arg(""), unitig_lengths_file_given(false),
    kunitigs_translation_file_arg(""), kunitigs_translation_file_given(false),
    overlaps_file_arg(""), overlaps_file_given(false),
    num_kunitigs_file_arg(""), num_kunitigs_file_given(false),
    max_nodes_allowed_arg(4000), max_nodes_allowed_given(false),
    output_arg("super_reads_output"), output_given(false),
    threads_arg(1), threads_given(false)
  { }

  joinKUnitigs_v3(int argc, char* argv[]) :
    min_overlap_length_arg(), min_overlap_length_given(false),
    mean_and_stdev_by_prefix_file_arg(""), mean_and_stdev_by_prefix_file_given(false),
    unitig_lengths_file_arg(""), unitig_lengths_file_given(false),
    kunitigs_translation_file_arg(""), kunitigs_translation_file_given(false),
    overlaps_file_arg(""), overlaps_file_given(false),
    num_kunitigs_file_arg(""), num_kunitigs_file_given(false),
    max_nodes_allowed_arg(4000), max_nodes_allowed_given(false),
    output_arg("super_reads_output"), output_given(false),
    threads_arg(1), threads_given(false)
  { parse(argc, argv); }

  void parse(int argc, char* argv[]) {
    static struct option long_options[] = {
      {"min-overlap-length", 1, 0, MIN_OVERLAP_LENGTH_OPT},
      {"mean-and-stdev-by-prefix-file", 1, 0, 'm'},
      {"unitig-lengths-file", 1, 0, 'u'},
      {"kunitigs-translation-file", 1, 0, KUNITIGS_TRANSLATION_FILE_OPT},
      {"overlaps-file", 1, 0, 'v'},
      {"num-kunitigs-file", 1, 0, NUM_KUNITIGS_FILE_OPT},
      {"max-nodes-allowed", 1, 0, MAX_NODES_ALLOWED_OPT},
      {"output", 1, 0, 'o'},
      {"threads", 1, 0, 't'},
      {"help", 0, 0, 'h'},
      {"usage", 0, 0, USAGE_OPT},
      {"version", 0, 0, 'V'},
      {0, 0, 0, 0}
    };
    static const char *short_options = "hVm:u:v:o:t:";

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
      case MIN_OVERLAP_LENGTH_OPT:
        min_overlap_length_given = true;
        min_overlap_length_arg = yaggo::conv_int<int>((const char*)optarg, err, false);
        CHECK_ERR(int_t, optarg, "    --min-overlap-length=int")
        break;
      case 'm':
        mean_and_stdev_by_prefix_file_given = true;
        mean_and_stdev_by_prefix_file_arg = optarg;
        break;
      case 'u':
        unitig_lengths_file_given = true;
        unitig_lengths_file_arg = optarg;
        break;
      case KUNITIGS_TRANSLATION_FILE_OPT:
        kunitigs_translation_file_given = true;
        kunitigs_translation_file_arg = optarg;
        break;
      case 'v':
        overlaps_file_given = true;
        overlaps_file_arg = optarg;
        break;
      case NUM_KUNITIGS_FILE_OPT:
        num_kunitigs_file_given = true;
        num_kunitigs_file_arg = optarg;
        break;
      case MAX_NODES_ALLOWED_OPT:
        max_nodes_allowed_given = true;
        max_nodes_allowed_arg = yaggo::conv_int<int>((const char*)optarg, err, false);
        CHECK_ERR(int_t, optarg, "    --max-nodes-allowed=int")
        break;
      case 'o':
        output_given = true;
        output_arg = optarg;
        break;
      case 't':
        threads_given = true;
        threads_arg = yaggo::conv_int<int>((const char*)optarg, err, false);
        CHECK_ERR(int_t, optarg, "-t, --threads=int")
        break;
      }
    }

    // Check that required switches are present
    if(!min_overlap_length_given)
      error("[    --min-overlap-length=int] required switch");
    if(!mean_and_stdev_by_prefix_file_given)
      error("[-m, --mean-and-stdev-by-prefix-file=path] required switch");
    if(!unitig_lengths_file_given)
      error("[-u, --unitig-lengths-file=path] required switch");
    if(!overlaps_file_given)
      error("[-v, --overlaps-file=path] required switch");
    if(!num_kunitigs_file_given)
      error("[    --num-kunitigs-file=path] required switch");

    // Parse arguments
    if(argc - optind != 1)
      error("Requires exactly 1 argument.");
    input_file_arg = argv[optind];
    ++optind;
  }

#define joinKUnitigs_v3_USAGE "Usage: joinKUnitigs_v3 [options] input-file:string"
  const char * usage() const { return joinKUnitigs_v3_USAGE; }
  void error(const char *msg) { 
    std::cerr << "Error: " << msg << "\n" << usage()
              << "\nUse --help for more information"
              << std::endl;
    exit(1);
  }

#define joinKUnitigs_v3_HELP "Join k-unitigs overlapping mate pairs of an insert.\n\nFor this exec we are using the unitig numbers starting from 0.\n\n" \
  "Options (default value in (), *required):\n" \
  "     --min-overlap-length=int            *Minimum length of an overlap between unitigs\n" \
  " -m, --mean-and-stdev-by-prefix-file=path\n                                         *File containing the mean and stdev for each prefix library.\n" \
  " -u, --unitig-lengths-file=path          *File containing the length of the unitigs.\n" \
  "     --kunitigs-translation-file=path     File containing map from original unitigs to new (longer) unitigs.\n" \
  " -v, --overlaps-file=path                *Celera-style overlap file between unitigs in binary format.\n" \
  "     --num-kunitigs-file=path            *File containing the number of k-unitigs.\n" \
  "     --max-nodes-allowed=int              Max records allowed when trying to join a mate pair. (4000)\n" \
  " -o, --output=string                      Output file (super_reads_output)\n" \
  " -t, --threads=int                        Number of threads (1)\n" \
  "     --usage                              Usage\n" \
  " -h, --help                               This message\n" \
  " -V, --version                            Version"

  const char * help() const { return joinKUnitigs_v3_HELP; }
#define joinKUnitigs_v3_HIDDEN "Hidden options:"

  const char * hidden() const { return joinKUnitigs_v3_HIDDEN; }
  void print_version(std::ostream &os = std::cout) const {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
    os << PACKAGE_VERSION << "\n";
  }
  void dump(std::ostream &os = std::cout) {
    os << "min_overlap_length_given:" << min_overlap_length_given << " min_overlap_length_arg:" << min_overlap_length_arg << "\n";
    os << "mean_and_stdev_by_prefix_file_given:" << mean_and_stdev_by_prefix_file_given << " mean_and_stdev_by_prefix_file_arg:" << mean_and_stdev_by_prefix_file_arg << "\n";
    os << "unitig_lengths_file_given:" << unitig_lengths_file_given << " unitig_lengths_file_arg:" << unitig_lengths_file_arg << "\n";
    os << "kunitigs_translation_file_given:" << kunitigs_translation_file_given << " kunitigs_translation_file_arg:" << kunitigs_translation_file_arg << "\n";
    os << "overlaps_file_given:" << overlaps_file_given << " overlaps_file_arg:" << overlaps_file_arg << "\n";
    os << "num_kunitigs_file_given:" << num_kunitigs_file_given << " num_kunitigs_file_arg:" << num_kunitigs_file_arg << "\n";
    os << "max_nodes_allowed_given:" << max_nodes_allowed_given << " max_nodes_allowed_arg:" << max_nodes_allowed_arg << "\n";
    os << "output_given:" << output_given << " output_arg:" << output_arg << "\n";
    os << "threads_given:" << threads_given << " threads_arg:" << threads_arg << "\n";
    os << "input_file_arg:" << input_file_arg << "\n";
  }
private:
};
#endif // __JOINKUNITIGS_V3_HPP__"
