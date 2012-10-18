#include <cstdio>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <charbuf.hpp>
#include <err.hpp>
#include <misc.hpp>
#include <src2/collectReadSequencesForLocalGapClosing_cmdline.hpp>

// THE NEXT 2 LINES ARE COPIED FROM GUILLAUME
// GLOBAL: command line switches
cmdline_parse args;

// Names of files and directories

struct numAndOriStruct {
     int groupNum;
     int ori;
};

typedef std::string stdString;
typedef std::vector<stdString> vectorOfStrings;
typedef std::set<stdString> setOfStrings;
typedef std::unordered_map<stdString, stdString> readNameToReadSequence;
typedef std::unordered_map<stdString, int> stringToIntMap;
typedef std::unordered_map<stdString, stdString> stringToStringMap;

stdString getReadMateName (stdString readName);
bool loadNeededReads (setOfStrings &readIsNeeded, readNameToReadSequence &readSeq);
void checkArgs (void);


int main(int argc, char **argv)
{
     FILE *infile; // , *outfile;
     std::ofstream outfile;
     charb line(100);
     vectorOfStrings fauxReadGroups;
     stringToIntMap group;
     std::unordered_map<stdString, char> end;
     setOfStrings readIsNeeded;
     readNameToReadSequence readSeq; // Read name to read sequence
     struct stat     statbuf;
     args.parse (argc, argv);
     checkArgs ();
     infile = fopen(args.faux_reads_file_arg, "r");
     while (fgets (line, 100, infile)) {
	  fauxReadGroups.push_back(stdString(line));
	  char *cptr = line+1;
	  char *fauxReadName;
	  char *saveptr;
	  fauxReadName = strtok_r(cptr, " \t\n", &saveptr);
	  stdString fauxReadNameStr = stdString(fauxReadName);
	  group[fauxReadNameStr] = fauxReadGroups.size()-1;
	  end[fauxReadNameStr] = 'F';
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += stdString (line);
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += stdString (line);
	  cptr = line+1;  // Don't know if it was reset (what happens in scanf)
	  fauxReadName = strtok_r(cptr, " \t\n", &saveptr);
	  fauxReadNameStr = stdString(fauxReadName);
	  group[fauxReadNameStr] = fauxReadGroups.size()-1;
	  end[fauxReadNameStr] = 'R';
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += stdString (line);
     }
     fclose (infile);

     int numFauxReadGroups = (int) fauxReadGroups.size();
     
     infile = fopen (args.faux_read_matches_to_kunis_file_arg, "r");
     stdString fauxRead;
     ExpBuffer<char *> flds;
     char matchToKeepMate;

     typedef std::vector<numAndOriStruct> numAndOriList;
     numAndOriList emptyNumAndOriList;
     std::vector<numAndOriList> kUniMatchStructs;
     std::vector<vectorOfStrings> readsInGroup, mateReadsInGroup;
     vectorOfStrings emptyStringVector;
     for (int i=0; i<numFauxReadGroups; i++) {
	  readsInGroup.push_back (emptyStringVector);
	  mateReadsInGroup.push_back (emptyStringVector);
     }

     int maxKUniSeen = -1;
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  // We will keep the info on what kind of a match we need to bring
	  // in a mate as well as a read
	  fauxRead = stdString (flds[0]);
	  int fauxReadGroup = group[fauxRead];
	  char fauxReadGapEnd = end[fauxRead];
	  for (unsigned int i=2; i<flds.size(); i+=3) {
	       int kUni = atoi (flds[i]);
	       char relOri = *flds[i+2];
	       while (kUni > maxKUniSeen) {
		    ++maxKUniSeen;
		    kUniMatchStructs.push_back(emptyNumAndOriList); }
	       if (fauxReadGapEnd == relOri)
		    matchToKeepMate = 'F';
	       else
		    matchToKeepMate = 'R';
	       numAndOriStruct numAndOri;
	       numAndOri.groupNum = fauxReadGroup;
	       numAndOri.ori = matchToKeepMate;
	       kUniMatchStructs[kUni].push_back(numAndOri);
	  }
     }
     fclose (infile);

     infile = fopen (args.read_matches_to_kunis_file_arg, "r");
     // Ex. line: pe2836 101 610 10 F
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  stdString readName = stdString (flds[0]);
	  int flds0len = strlen(flds[0]);
	  if (flds[0][flds0len-1] % 2 == 0)
	       ++flds[0][flds0len-1];
	  else
	       --flds[0][flds0len-1];
	  stdString readMate = stdString (flds[0]);
	  int kUni = atoi (flds[2]);
	  char relOri = *(flds[4]);
	  for (unsigned int i=0; i<kUniMatchStructs[kUni].size(); i++) {
	       readsInGroup[kUniMatchStructs[kUni][i].groupNum].push_back(readName);
	       readIsNeeded.insert (readName);
	       if (kUniMatchStructs[kUni][i].ori == relOri) {
		    readIsNeeded.insert (readMate);
		    mateReadsInGroup[kUniMatchStructs[kUni][i].groupNum].push_back(readMate); }
	  }
     }

     if (stat (args.dir_for_gaps_arg, &statbuf) != 0) {
	  charb cmd(100);
	  sprintf (cmd, "mkdir %s", args.dir_for_gaps_arg);
	  fprintf (stderr, "%s\n", (char *) cmd);
	  system (cmd);
     }
     int outputGroupNum = 0;
     while (1) {
//	  bool isAtEnd = loadNeededReads();
	  bool isAtEnd = loadNeededReads (readIsNeeded, readSeq);
	  if (readSeq.empty())
	       break;
	  charb outfileName(100);
	  sprintf (outfileName, "%s/readFile.%03d", args.dir_for_gaps_arg, outputGroupNum);
	  outfile.open (outfileName);
	  if (args.output_dir_progress_flag)
	       fprintf (stderr, "Outputting file %s\n", (char *) outfileName);
	  
	  for (unsigned int grp=0; grp<readsInGroup.size(); grp++) {
	       charb istr2(20), newDir(20);
	       outfile << '>' << grp << "\n";
	       
	       stringToIntMap willBeOutput;
	       willBeOutput.clear();
	       vectorOfStrings readsToOutput;
	       readsToOutput.clear();
	       for (unsigned int readNum=0; readNum<readsInGroup[grp].size(); readNum++) {
		    stdString readName = readsInGroup[grp][readNum];
		    stringToIntMap::iterator it = willBeOutput.find (readName);
		    if (it != willBeOutput.end())
			 continue;
		    readsToOutput.push_back (readName);
		    willBeOutput[readName] = 2; }
	       for (unsigned int readNum=0; readNum<mateReadsInGroup[grp].size(); readNum++) {
		    stdString readName = mateReadsInGroup[grp][readNum];
		    stringToIntMap::iterator it = willBeOutput.find (readName);
		    if (it != willBeOutput.end())
			 continue;
		    willBeOutput[readName] = 1; }
	       setOfStrings alreadyOutput;
	       alreadyOutput.clear();
	       vectorOfStrings bothMatesHaveMatches, mateBroughtInViaMatePair, readOnlys;
	       bothMatesHaveMatches.clear();
	       mateBroughtInViaMatePair.clear();
	       readOnlys.clear();
	       for (unsigned int readNum=0; readNum<readsToOutput.size(); readNum++) {
		    stdString readName = readsToOutput[readNum];
		    if (alreadyOutput.find(readName) != alreadyOutput.end())
			 continue;
		    stdString readMate = getReadMateName (readName);
		    alreadyOutput.insert (readName);
		    stringToIntMap::iterator it2 = willBeOutput.find (readMate);
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
		    stringToStringMap::iterator readSeqIt =
			 readSeq.find (bothMatesHaveMatches[readNum]);
		    if (readSeqIt == readSeq.end())
			 continue;
		    outfile << bothMatesHaveMatches[readNum] << " B\n" <<
			 readSeqIt->second << "\n";
	       }
	       int outCount = 0;
	       for (unsigned int readNum=0; readNum<mateBroughtInViaMatePair.size(); readNum++) {
		    ++outCount;
		    stdString extraStr;
		    stdString readName = mateBroughtInViaMatePair[readNum];
		    stringToStringMap::iterator readSeqIt =
			 readSeq.find (readName);
		    if (readSeqIt == readSeq.end())
			 continue;
		    if (outCount % 2 == 1)
			 extraStr = stdString ("match");
		    else
			 extraStr = stdString ("mate");
		    outfile << readName << ' ' << extraStr << " M\n" <<
			 readSeq[readName] << "\n";
	       }
	       for (unsigned int readNum=0; readNum<readOnlys.size(); readNum++) {
		    stdString readName = readOnlys[readNum];
		    stringToStringMap::iterator readSeqIt =
			 readSeq.find (readName);
		    if (readSeqIt == readSeq.end())
			 continue;
		    stdString mateRead = getReadMateName (readName);
		    outfile << readName << " O\n" << readSeq[readName] << '\n'
			    << mateRead << " N\n" << "N\n";
	       }
	  }
	  outfile.close();
	  ++outputGroupNum;
	       
	  if (isAtEnd)
	       break;
     }

     return (0);
}

bool loadNeededReads (setOfStrings &readIsNeeded, readNameToReadSequence &readSeq)
{
     FILE *infile;
     static unsigned int currentFileNum = 0;
     static off_t currentFileOffset = (off_t) 0;
     const char *localReadsFile;
     uint64_t numReadSeqsLoaded;
     
     readSeq.clear();
     numReadSeqsLoaded = 0;
     while (currentFileNum < args.reads_file_arg.size()) {
	  localReadsFile = args.reads_file_arg[currentFileNum];
	  infile = fopen (localReadsFile, "r");
	  fseeko (infile, currentFileOffset, SEEK_SET);
	  int offsetStart = strlen (localReadsFile)-5;
	  charb line(100), readNameStr(100);
	  char *cptr = line+1;
	  stdString readName;
	  bool isNeeded = false;
	  if (offsetStart < 0) offsetStart = 0;
	  const char *cptr2 = localReadsFile + offsetStart;
	  if (strcmp (cptr2, "fastq") == 0)
	       goto fastQ;
	  while (fgets (line, 100, infile)) {
	       int lineLen;
	       cptr = line+1;
	       if (line[0] == '>') {
		    char *saveptr;
		    lineLen = strlen(line);
		    readNameStr = strtok_r(cptr, " \t\n", &saveptr);
		    readName = stdString (readNameStr);
		    if (readIsNeeded.find(readName) != readIsNeeded.end())
			 isNeeded = true;
		    else
			 isNeeded = false;
		    // What to do if we should stop before this read
		    if ((((numReadSeqsLoaded == args.max_reads_in_memory_arg-1) && (readNameStr[strlen(readNameStr)-1] % 2 == 0)) ||
			 (numReadSeqsLoaded == args.max_reads_in_memory_arg)) && isNeeded) {
			 currentFileOffset = ftello (infile);
			 currentFileOffset -= lineLen;
			 fclose (infile);
			 return false;
		    }
		    continue; }
	       if (isNeeded) {
		    line[strlen(line)-1] = 0; // chop (line);
		    if (readSeq.find(readName) == readSeq.end()) {
			 readSeq[readName] = stdString ("");
			 ++numReadSeqsLoaded; }
		    readSeq[readName] += stdString (line);
	       }
	  }
	  goto endProcessingReads;
     fastQ:
	  while (fgets (line, 100, infile)) {
	       int lineLen;
	       lineLen = strlen (line);
	       cptr = line+1;
	       char *saveptr;
	       isNeeded = false;
	       readNameStr = strtok_r(cptr, " \t\n", &saveptr);
	       readName = stdString (readNameStr);
	       if (readIsNeeded.find(readName) != readIsNeeded.end())
		    isNeeded = true;
	       else
		    isNeeded = false;
	       // What to do if we should stop before this read
	       if ((((numReadSeqsLoaded == args.max_reads_in_memory_arg-1) && (readNameStr[strlen(readNameStr)-1] % 2 == 0)) ||
		    (numReadSeqsLoaded == args.max_reads_in_memory_arg)) && isNeeded) {
		    currentFileOffset = ftello (infile);
		    currentFileOffset -= lineLen;
		    fclose (infile);
		    return false;
	       }

	       fgets (line, 100, infile);
	       line[strlen(line)-1] = 0; // chop (line);
	       if (isNeeded) {
		    readSeq[readName] = stdString (line);
		    ++numReadSeqsLoaded;
	       }
	       fgets (line, 100, infile); // Skip quality bases
	       fgets (line, 100, infile);
	  }
     endProcessingReads:
	  fclose (infile);
	  ++currentFileNum;
	  currentFileOffset = (off_t) 0;
     }
     return true;
}

stdString getReadMateName (stdString readName)
{
     stdString readMate;
     int readNameLen = readName.length();
     charb readNameStr(100);
     if (readNameLen == 0)
	  return (stdString (""));
     strcpy (readNameStr, readName.c_str());
     if (readNameStr[readNameLen-1] % 2 == 0)
	       ++readNameStr[readNameLen-1];
	  else
	       --readNameStr[readNameLen-1];
     readMate = stdString (readNameStr);

     return (readMate);
}

void checkArgs (void)
{
     struct stat     statbuf;
     int fail = 0;
     if (stat (args.faux_reads_file_arg, &statbuf) != 0) {
	  std::cerr << "Faux reads file '" << args.faux_reads_file_arg << "' doesn't exist.\n";
	  fail = 1; }
     if (stat (args.faux_read_matches_to_kunis_file_arg, &statbuf) != 0) {
	  std::cerr << "File of faux read matches to k-unitigs '" << args.faux_read_matches_to_kunis_file_arg << "' doesn't exist.\n";
	  fail = 1; }
     if (stat (args.read_matches_to_kunis_file_arg, &statbuf) != 0) {
	  std::cerr << "File of read matches to k-unitigs '" << args.read_matches_to_kunis_file_arg << "' doesn't exist.\n";
	  fail = 1; }
     for (unsigned int i=0; i<args.reads_file_arg.size(); i++) {
	  const char *readFile = args.reads_file_arg[i];
	  if (stat (readFile, &statbuf) != 0) {
	       std::cerr << "Read file '" << readFile << "' doesn't exist.\n";
	       fail = 1; }
     }
     
     if (fail)
	  exit (1);
}


     
