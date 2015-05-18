#include <cstdio>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <charbuf.hpp>
#include <jellyfish/err.hpp>
#include <misc.hpp>
#include <src2/collectReadSequencesForLocalGapClosing_cmdline.hpp>

// GLOBAL: command line switches
cmdline_parse args;

struct numAndOriStruct {
     int groupNum;
     int ori;
};

struct group_info {
  int  group;
  char end;
};

// Match information is kept for each mate pairs, indexed by the name
// of the first read in a pair. For each pair we keep for each group
// the reason why the first read in a pair or its mate might be
// needed: it has a match, the mate has a match, it is pulled in by
// its mate, the mate is pulled in by the read. If any of those flag
// is true in any of the group, then the needed flag is ORed
// accordingly to isNeeded and mateNeeded.
enum needed_type : char { noneNeeded = 0, isNeeded = 1, mateNeeded = 2 };
enum match_type : char { noMatch = 0, hasMatch = 1, mateMatch = 2, pulledIn = 4, matePulled = 8 };
struct group_matches_type {
  int  groupNum;
  char match;
};
struct pair_group_list_type {
  char                            needed;
  std::vector<group_matches_type> matches;
};
typedef std::unordered_map<std::string, pair_group_list_type> pair_group_info_type;


typedef std::string stdString;
typedef std::vector<stdString> vectorOfStrings;
typedef std::unordered_set<stdString> setOfStrings;
typedef std::unordered_map<stdString, stdString> readNameToReadSequence;
typedef std::unordered_map<stdString, int> stringToIntMap;
typedef std::unordered_map<stdString, stdString> stringToStringMap;
typedef std::vector<numAndOriStruct> numAndOriList;

struct header_sequence_type {
  std::string header;
  std::string sequence;
};
struct read_list_type {
  size_t index; // index in readSeq
  char   match; // match type
};


inline bool readIsFirstOfPair (const stdString& readName);
stdString getReadMateName (const stdString& readName);
bool loadNeededReads (const pair_group_info_type& pair_info, std::vector<header_sequence_type>& readSeq,
                      std::vector<std::vector<read_list_type> >& reads_in_group);
void checkArgs (void);
std::unordered_map<stdString, group_info> load_group_info(FILE* infile);
std::vector<numAndOriList> load_read_matches_unis(FILE* infile, const std::unordered_map<stdString, group_info>& groups);
pair_group_info_type load_read_mate_groups(FILE* infile, const std::vector<numAndOriList>& kUniMatchStructs);
void read_and_write_reads(const pair_group_info_type& pair_info, const unsigned int numGrps);


int main(int argc, char **argv)
{
     args.parse (argc, argv);
     const unsigned int numGrpsInGrp = args.num_joins_per_directory_arg;

     setOfStrings readIsNeeded;
     pair_group_info_type pair_info;
     unsigned int numGrps = 0;
     {
         FILE *infile;
         infile = fopen(args.faux_reads_file_arg, "r");
         if(!infile)
             cmdline_parse::error() << "Failed to open '" << args.faux_reads_file_arg << "' for reading";
         const std::unordered_map<stdString, group_info> groups = load_group_info(infile);
         fclose (infile);
         const unsigned int numFauxReadGroups = groups.size() / 2;
         numGrps = numFauxReadGroups / numGrpsInGrp + (numFauxReadGroups % numGrpsInGrp != 0);

         infile = fopen(args.faux_read_matches_to_kunis_file_arg, "r");
         if(!infile)
             cmdline_parse::error() << "Failed top open '" << args.faux_read_matches_to_kunis_file_arg << "' for reading";
         std::vector<numAndOriList> kUniMatchStructs = load_read_matches_unis(infile, groups);
         fclose(infile);

         infile = fopen (args.read_matches_to_kunis_file_arg, "r");
         if(!infile)
             cmdline_parse::error() << "Feiled to open '" << args.read_matches_to_kunis_file_arg << "' for reading";
         pair_info = load_read_mate_groups(infile, kUniMatchStructs);
         fclose(infile);
     }

     int mkdir_err = mkdir(args.dir_for_gaps_arg,
                           S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);
     if(mkdir_err != 0 && errno != EEXIST)
       cmdline_parse::error() << "Failed to create directory for gaps '" << args.dir_for_gaps_arg << "'";

     read_and_write_reads(pair_info, numGrps);

     return EXIT_SUCCESS;
}

inline void ignore_line(std::istream& is) {
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

bool loadNeededReads (const pair_group_info_type& pair_info, std::vector<header_sequence_type>& readSeq,
                      std::vector<std::vector<read_list_type> >& reads_in_group) {
    static unsigned int  currentFileNum    = 0;
    static bool          currentFileFasta  = false;
    static std::ifstream infile;

    for(auto& group : reads_in_group)
      group.clear();

    charb  seq_tmp;
    size_t read_index = 0;
    bool   done       = false;

    while (!done) {
      if(!infile.is_open() || infile.peek() == EOF) {
        infile.close();
        if(currentFileNum >= args.reads_file_arg.size())
          break;
        infile.open(args.reads_file_arg[currentFileNum++]);
        switch(infile.peek()) {
        case '>': currentFileFasta = true; break;
        case '@': currentFileFasta = false; break;
        default: cmdline_parse::error() << "Invalid format. Expected '>' or '@'. Got '" << (char)infile.peek() << "'";
        }
      }

      for(int c = infile.get(); c != EOF; c = infile.get()) {
        if(c != (currentFileFasta ? '>' : '@'))
          cmdline_parse::error() << "Invalid file format. Expected '" << (currentFileFasta ? '>' : '@')
                                 << "' got '" << (char)c << "'";

        std::string& header = readSeq[read_index].header;
        getline(infile, header);
        const size_t line_size = header.size();
        const auto space = header.find_first_of(" \t\n");
        if(space != std::string::npos) header.resize(space);
        // Got enough reads?
        if ((read_index == args.max_reads_in_memory_arg - 1 && readIsFirstOfPair(header)) ||
            (read_index == args.max_reads_in_memory_arg)) {
          infile.seekg(-(line_size + 2), std::ios::cur); // 2 == \n + ('>' | '@')
          done = true;
          break;
        }

        const bool        isFirst   = readIsFirstOfPair(header);
        const std::string firstName = isFirst ? header : getReadMateName(header);
        const auto it = pair_info.find(firstName);
        if(it == pair_info.cend() || ((it->second.needed & (isFirst ? isNeeded : mateNeeded)) == 0)) {
          // File is ignored, skip lines
          if(currentFileFasta) {
            for(c = infile.peek(); c != '>' && c != EOF; c = infile.peek())
              ignore_line(infile);
          } else {
            for(int i = 0; i < 3; ++i)
              ignore_line(infile);
          }
        } else {
          // Record read in readSeq
          std::string& sequence = readSeq[read_index].sequence;
          if(currentFileFasta) { // potential multiline sequence
            getline(infile, seq_tmp);
            for(c = infile.peek(); c != '>' && c != EOF; c = infile.peek())
              getline_append(infile, seq_tmp);
            sequence.assign(seq_tmp, seq_tmp.size());
          } else {
            getline(infile, sequence);
            ignore_line(infile);
            ignore_line(infile);  // Quality scores
          }
          for(const auto& match : it->second.matches) {
            if((match.match & (isFirst ? (hasMatch | pulledIn) : (mateMatch | matePulled))) != 0)
              reads_in_group[match.groupNum].push_back({ read_index, match.match });
          }
          ++read_index;
        }
      }
    }
    return read_index > 0;
}

bool readIsFirstOfPair (const stdString& readName) {
  return readName.length() > 0 && (readName.back() % 2 == 0);
}

stdString getReadMateName (const stdString& readName) {
  stdString readMate(readName);
  if(readMate.length() > 0)
    readMate.back() += 1 - 2 * (readMate.back() % 2);

  return readMate;
}

std::unordered_map<stdString, group_info> load_group_info(FILE* infile) {
  std::unordered_map<stdString, group_info> groups;
  charb                                     line;

  for(int group_index = 0 ; fgets(line, infile); ++group_index) {
    char        *saveptr;
    const char*  fauxReadName     = strtok_r(line + 1, " \t\n", &saveptr);
    groups[fauxReadName] = { group_index, 'F' };
    fgets (line, infile);
    fgets (line, infile);
    fauxReadName                  = strtok_r(line + 1, " \t\n", &saveptr);
    groups[fauxReadName] = { group_index, 'R' };
    fgets (line, infile);
  }
  return groups;
}

std::vector<numAndOriList> load_read_matches_unis(FILE* infile, const std::unordered_map<stdString, group_info>& groups) {
  std::vector<numAndOriList> kUniMatchStructs;
  ExpBuffer<char *>          flds;
  charb                      line;

  while (fgets (line, infile)) {
    getFldsFromLine (line, flds);
    // We will keep the info on what kind of a match we need to bring
    // in a mate as well as a read
    const auto fauxReadInfo = groups.find(flds[0]);
    if(fauxReadInfo == groups.cend())
      cmdline_parse::error() << "Unknown group '" << flds[0] << "' in faux read matches file";
    for (unsigned int i=2; i<flds.size(); i+=3) {
      const int  kUni   = atoi (flds[i]);
      const char relOri = *flds[i+2];
      if(kUni >= (int)kUniMatchStructs.size())
        kUniMatchStructs.resize(kUni + 1);
      const char matchToKeepMate = fauxReadInfo->second.end == relOri ? 'F' : 'R';
      kUniMatchStructs[kUni].push_back({ fauxReadInfo->second.group, matchToKeepMate });
    }
  }
  return kUniMatchStructs;
}

pair_group_info_type load_read_mate_groups(FILE* infile, const std::vector<numAndOriList>& kUniMatchStructs) {
  pair_group_info_type res;
  charb                line;
  ExpBuffer<char *>    flds;

  // Ex. line: pe2836 101 610 10 F
  while(fgets (line, infile)) {
    getFldsFromLine (line, flds);
    const stdString& readName  = flds[0];
    const int        kUni      = atoi (flds[2]);
    const char       relOri    = *(flds[4]);
    const bool       isFirst   = readIsFirstOfPair(readName);
    const stdString  firstName = isFirst ? readName : getReadMateName(readName);

    auto& pair_info = res[firstName];
    const char my_match     = isFirst ? hasMatch : mateMatch;
    const char my_pulled    = isFirst ? matePulled : pulledIn;
    bool       other_needed = false;
    if(pair_info.needed == noneNeeded) { // Means first time see this mate pair -> copy match information
      pair_info.matches.reserve(kUniMatchStructs[kUni].size());
      for(const auto& match_info : kUniMatchStructs[kUni]) {
        const bool pulled  = match_info.ori == relOri;
        other_needed      |= pulled;
        pair_info.matches.push_back({ match_info.groupNum / args.num_joins_per_directory_arg,
              (char)(my_match | (pulled ? my_pulled : noMatch)) });
      }
    } else { // Already seen the pair -> merge match information
      std::unordered_map<int, char> group2match;
      for(const auto& match_info : kUniMatchStructs[kUni]) { // put in hash for fast access
        const bool pulled  = match_info.ori == relOri;
        other_needed      |= pulled;
        group2match[match_info.groupNum / args.num_joins_per_directory_arg] = my_match | (pulled ? my_pulled : noMatch);
      }
      for(auto& group_match : pair_info.matches) { // do merge
        auto it = group2match.find(group_match.groupNum);
        if(it != group2match.end()) {
          group_match.match |= it->second;
          group2match.erase(it);
        }
      }
      pair_info.matches.reserve(pair_info.matches.size() + group2match.size());
      for(const auto& match_info : group2match) { // add remaining
        pair_info.matches.push_back({ match_info.first, match_info.second });
      }
    }
    pair_info.needed |= (isFirst ? isNeeded : mateNeeded);
    if(other_needed)
      pair_info.needed |= isFirst ? mateNeeded : isNeeded;

  }
  return res;
}

inline char match_str(char x) {
  if((x & (hasMatch | mateMatch)) == (hasMatch | mateMatch))
     return 'B';
  if((x & (pulledIn | matePulled)) != 0)
    return 'M';
  return 'O';
}

void read_and_write_reads(const pair_group_info_type& pair_info, const unsigned int numGrps) {
  std::vector<header_sequence_type> readSeq(args.max_reads_in_memory_arg + 1);
  std::vector<std::vector<read_list_type> > reads_in_group(numGrps);
  charb outfileName;

  for(unsigned int outputGroupNum = 0; loadNeededReads(pair_info, readSeq, reads_in_group); ++outputGroupNum) {
    sprintf (outfileName, "%s/readFile.%03d", args.dir_for_gaps_arg, outputGroupNum);
    std::ofstream outfile(outfileName);

    for(unsigned int grp = 0; grp < numGrps; ++grp) {
      outfile << '>' << grp << '\n';
      for(const auto& read : reads_in_group[grp])
        outfile << readSeq[read.index].header << ' ' << grp
                << ' ' << match_str(read.match) << '\n'
                << readSeq[read.index].sequence << '\n';
    }
  }
}
