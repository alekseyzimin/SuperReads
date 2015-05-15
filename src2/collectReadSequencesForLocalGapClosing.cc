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

typedef std::string stdString;
typedef std::vector<stdString> vectorOfStrings;
typedef std::unordered_set<stdString> setOfStrings;
typedef std::unordered_map<stdString, stdString> readNameToReadSequence;
typedef std::unordered_map<stdString, int> stringToIntMap;
typedef std::unordered_map<stdString, stdString> stringToStringMap;
typedef std::vector<numAndOriStruct> numAndOriList;

inline bool readIsFirstOfPair (const stdString& readName);
stdString getReadMateName (const stdString& readName);
void loadNeededReads (const setOfStrings &readIsNeeded, readNameToReadSequence &readSeq);
void checkArgs (void);
std::unordered_map<stdString, group_info> load_group_info(FILE* infile);
std::vector<numAndOriList> load_read_matches_unis(FILE* infile, const std::unordered_map<stdString, group_info>& groups);
setOfStrings load_read_mate_groups(FILE* infile, const std::vector<numAndOriList>& kUniMatchStruct,
                                   std::vector<vectorOfStrings>& readsInGroup, std::vector<vectorOfStrings>& mateReadsInGroup);
void read_and_write_reads(const setOfStrings& readIsNeeded,
                          const std::vector<vectorOfStrings>& readsInGroup,
                          const std::vector<vectorOfStrings>& mateReadsInGroup);


int main(int argc, char **argv)
{
     args.parse (argc, argv);

     setOfStrings readIsNeeded;
     std::vector<vectorOfStrings> readsInGroup, mateReadsInGroup;
     {
         FILE *infile;
         infile = fopen(args.faux_reads_file_arg, "r");
         if(!infile)
             cmdline_parse::error() << "Failed to open '" << args.faux_reads_file_arg << "' for reading";
         const std::unordered_map<stdString, group_info> groups = load_group_info(infile);
         fclose (infile);
         int numFauxReadGroups = (int)groups.size() / 2;

         infile = fopen(args.faux_read_matches_to_kunis_file_arg, "r");
         if(!infile)
             cmdline_parse::error() << "Failed top open '" << args.faux_read_matches_to_kunis_file_arg << "' for reading";
         std::vector<numAndOriList> kUniMatchStructs = load_read_matches_unis(infile, groups);
         fclose(infile);

         infile = fopen (args.read_matches_to_kunis_file_arg, "r");
         if(!infile)
             cmdline_parse::error() << "Feiled to open '" << args.read_matches_to_kunis_file_arg << "' for reading";
         readsInGroup.resize(numFauxReadGroups);
         mateReadsInGroup.resize(numFauxReadGroups);
         readIsNeeded = load_read_mate_groups(infile, kUniMatchStructs, readsInGroup, mateReadsInGroup);
         fclose(infile);
     }

     int mkdir_err = mkdir(args.dir_for_gaps_arg,
                           S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);
     if(mkdir_err != 0 && errno != EEXIST)
       cmdline_parse::error() << "Failed to create directory for gaps '" << args.dir_for_gaps_arg << "'";

     read_and_write_reads(readIsNeeded, readsInGroup, mateReadsInGroup);

     return EXIT_SUCCESS;
}

inline void ignore_line(std::istream& is) {
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void loadNeededReads (const setOfStrings &readIsNeeded, readNameToReadSequence &readSeq)
{
    static unsigned int  currentFileNum    = 0;
    static off_t         currentFileOffset = (off_t) 0;
    static bool          currentFileFasta  = false;

    readSeq.clear();
    while (currentFileNum < args.reads_file_arg.size()) {
        std::ifstream infile(args.reads_file_arg[currentFileNum]);
        if(currentFileOffset > 0) {
            infile.seekg(currentFileOffset);
        } else {
            switch(infile.peek()) {
            case '>': currentFileFasta = true; break;
            case '@': currentFileFasta = false; break;
            default: cmdline_parse::error() << "Invalid format. Expected '>' or '@'. Got '" << (char)infile.peek() << "'";
            }
        }
        charb        header, sequence;
        const char*  readNameStr = nullptr;
        char        *saveptr;
        for(int c = infile.peek(); c != EOF; c = infile.peek()) {
            getline(infile, header);
            readNameStr = strtok_r(header + 1, " \t\n", &saveptr);
            if(c != (currentFileFasta ? '>' : '@'))
                cmdline_parse::error() << "Invalid file format. Expected '" << (currentFileFasta ? '>' : '@')
                                       << "' got '" << (char)c << "'";
            // Got enough reads?
            if ((readSeq.size() == args.max_reads_in_memory_arg-1 && readNameStr[strlen(readNameStr) - 1] % 2 == 0) ||
                (readSeq.size() == args.max_reads_in_memory_arg)) {
                currentFileOffset = (off_t)infile.tellg() - header.size() - 1;
                return;
            }
            getline(infile, sequence);
            if(currentFileFasta) { // potential multiline sequence
                for(c = infile.peek(); c != '>' && c != EOF; c = infile.peek())
                    getline_append(infile, sequence);
            }
            if(readIsNeeded.find(readNameStr) != readIsNeeded.end())
                readSeq[readNameStr] = sequence;
            if(!currentFileFasta) { // Quality scores
                ignore_line(infile);
                ignore_line(infile);
            }
        }
        currentFileOffset = (off_t) 0;
        ++currentFileNum;
    }
}

bool readIsFirstOfPair (const stdString& readName)
{
  return readName.length() > 0 && (readName.back() % 2 == 0);
}

stdString getReadMateName (const stdString& readName)
{
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

setOfStrings load_read_mate_groups(FILE* infile, const std::vector<numAndOriList>& kUniMatchStructs,
                                   std::vector<vectorOfStrings>& readsInGroup, std::vector<vectorOfStrings>& mateReadsInGroup) {
  setOfStrings readIsNeeded;
  charb line;
  ExpBuffer<char *> flds;

  // Ex. line: pe2836 101 610 10 F
  while (fgets (line, infile)) {
    getFldsFromLine (line, flds);
    const stdString readName = flds[0];
    const stdString readMate = getReadMateName(readName);
    const int       kUni     = atoi (flds[2]);
    const char      relOri   = *(flds[4]);
    for (unsigned int i=0; i<kUniMatchStructs[kUni].size(); i++) {
      readsInGroup[kUniMatchStructs[kUni][i].groupNum].push_back(readName);
      readIsNeeded.insert (readName);
      if (kUniMatchStructs[kUni][i].ori == relOri) {
        readIsNeeded.insert (readMate);
        mateReadsInGroup[kUniMatchStructs[kUni][i].groupNum].push_back(readMate); }
    }
  }
  return readIsNeeded;
}

void read_and_write_reads(const setOfStrings& readIsNeeded,
                          const std::vector<vectorOfStrings>& readsInGroup,
                          const std::vector<vectorOfStrings>& mateReadsInGroup) {
    int                      outputGroupNum = 0;
    std::map<stdString, int> readType;
    readNameToReadSequence   readSeq; // Read name to read sequence
    std::unordered_map<stdString, stdString>  outputReadHdr;
    while (true) {
        loadNeededReads (readIsNeeded, readSeq);
        if (readSeq.empty())
            break;
        charb outfileName(100);
        sprintf (outfileName, "%s/readFile.%03d", args.dir_for_gaps_arg, outputGroupNum);
        std::ofstream outfile(outfileName);

        const int numGrpsInGrp = args.num_joins_per_directory_arg;
        for (unsigned int grp=0; grp<readsInGroup.size(); grp++) {
            stringToIntMap willBeOutput;
            vectorOfStrings readsToOutput;
            for(const auto& readName : readsInGroup[grp]) {
                if (willBeOutput.find (readName) != willBeOutput.end())
                    continue;
                readsToOutput.push_back (readName);
                willBeOutput[readName] = 2; }
            for (const auto& readName : mateReadsInGroup[grp]) {
                if (willBeOutput.find (readName) != willBeOutput.end())
                    continue;
                willBeOutput[readName] = 1; }
            setOfStrings alreadyOutput;
            vectorOfStrings bothMatesHaveMatches, mateBroughtInViaMatePair, readOnlys;
            for(const auto& readName : readsToOutput) {
                if (alreadyOutput.find(readName) != alreadyOutput.end())
                    continue;
                const stdString readMate = getReadMateName (readName);
                alreadyOutput.insert (readName);
                auto it2 = willBeOutput.find (readMate);
                if (it2 != willBeOutput.end()) {
                    alreadyOutput.insert(readMate);
                    if (it2->second == 2) {
                        bothMatesHaveMatches.push_back (readName);
                        bothMatesHaveMatches.push_back (readMate); }
                    else {
                        mateBroughtInViaMatePair.push_back (readName);
                        mateBroughtInViaMatePair.push_back (readMate); }
                }
                else
                    readOnlys.push_back (readName);
            }
            for (unsigned int readNum=0; readNum<bothMatesHaveMatches.size(); readNum++) {
                if (readSeq.find (bothMatesHaveMatches[readNum]) == readSeq.end())
                    continue;
                outputReadHdr[bothMatesHaveMatches[readNum]] += (" " + std::to_string(grp) + " B");
                readType[bothMatesHaveMatches[readNum]] = 1;
            }
            int outCount = 0;
            for (const auto& readName : mateBroughtInViaMatePair) {
                ++outCount;
                if (readSeq.find (readName) == readSeq.end())
                    continue;
                readType[readName] = 1;
                outputReadHdr[readName] += (" " + std::to_string(grp) + " M");
            }
            for (const auto& readName : readOnlys) {
                if (readSeq.find (readName) == readSeq.end())
                    continue;
                const stdString mateRead = getReadMateName (readName);
                if (readType.find (readName) == readType.end())
                    readType[readName] = 2;
                else if (readType[readName] > 2)
                    readType[readName] = 2;
                outputReadHdr[readName] += (" " + std::to_string(grp) + " O");
                if (readType.find (mateRead) == readType.end())
                    readType[mateRead] = 3;
                outputReadHdr[mateRead] += (" " + std::to_string(grp) + " N");
            }
            if ((grp+1) % numGrpsInGrp == 0 || (grp+1) == readsInGroup.size()) {
                outfile << '>' << (grp / numGrpsInGrp) << "\n";
                for (const auto& readInfo : readType) {
                    const stdString& readName = readInfo.first;
                    if (! readIsFirstOfPair (readName))
                        continue;
                    const stdString mateRead = getReadMateName (readName);
                    if (readType[readName] != 3)
                        outfile << readName << outputReadHdr[readName] << '\n' << readSeq[readName] << '\n';
                    if (readType[mateRead] != 3)
                        outfile << mateRead << outputReadHdr[mateRead] << '\n' << readSeq[mateRead] << '\n';
                }
                outputReadHdr.clear();
                readType.clear();
            }
        }
        outfile.close();
        ++outputGroupNum;
    }
}
