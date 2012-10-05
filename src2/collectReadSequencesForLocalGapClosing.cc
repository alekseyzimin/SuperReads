#include <cstdio>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <charbuf.hpp>
#include <err.hpp>
#include <misc.hpp>

// Names of files and directories
std::string fauxReadsFile, fauxReadMatchesToKUnisFile, readMatchesToKUnisFile, dirForGaps;
charb bothHaveMatchFile(100), mateBroughtInViaMatePairFile(100),
     readsOnlyFile(100), newReadsFile(100);
std::set<std::string> readIsNeeded;
std::map<std::string, std::string> readSeq; // Read name to read sequence
bool outputDirProgress;

struct numAndOriStruct {
     int groupNum;
     int ori;
};
std::vector<std::string> readFiles;
int64_t maxReadsInMemory;

std::string getReadMateName (std::string readName);
bool loadNeededReads (void);
void processArgs (int argc, char **argv);

int main(int argc, char **argv)
{
     FILE *infile; // , *outfile;
     std::ofstream outfile;
     charb line(100);
     std::vector<std::string> fauxReadGroups;
     std::map<std::string,int> group;
     std::map<std::string,char> end;
     struct stat     statbuf;
     processArgs (argc, argv);
     infile = fopen(fauxReadsFile.c_str(), "r");
     while (fgets (line, 100, infile)) {
	  fauxReadGroups.push_back(std::string(line));
	  char *cptr = line+1;
	  char *fauxReadName;
	  char *saveptr;
	  fauxReadName = strtok_r(cptr, " \t\n", &saveptr);
	  std::string fauxReadNameStr = std::string(fauxReadName);
	  group[fauxReadNameStr] = fauxReadGroups.size()-1;
	  end[fauxReadNameStr] = 'F';
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += std::string (line);
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += std::string (line);
	  cptr = line+1;  // Don't know if it was reset (what happens in scanf)
	  fauxReadName = strtok_r(cptr, " \t\n", &saveptr);
	  fauxReadNameStr = std::string(fauxReadName);
	  group[fauxReadNameStr] = fauxReadGroups.size()-1;
	  end[fauxReadNameStr] = 'R';
	  fgets (line, 100, infile);
	  fauxReadGroups.back() += std::string (line);
     }
     fclose (infile);

     int numFauxReadGroups = (int) fauxReadGroups.size();
     
     infile = fopen (fauxReadMatchesToKUnisFile.c_str(), "r");
     std::string fauxRead;
     ExpBuffer<char *> flds;
     char matchToKeepMate;

     typedef std::vector<numAndOriStruct> numAndOriList;
     numAndOriList emptyNumAndOriList;
     std::vector<numAndOriList> kUniMatchStructs;
     std::vector<std::vector<std::string> > readsInGroup, mateReadsInGroup;
     std::vector<std::string> emptyStringVector;
     for (int i=0; i<numFauxReadGroups; i++) {
	  readsInGroup.push_back (emptyStringVector);
	  mateReadsInGroup.push_back (emptyStringVector);
     }

     int maxKUniSeen = -1;
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  // We will keep the info on what kind of a match we need to bring
	  // in a mate as well as a read
	  fauxRead = std::string (flds[0]);
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

     infile = fopen (readMatchesToKUnisFile.c_str(), "r");
     // Ex. line: pe2836 101 610 10 F
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  std::string readName = std::string (flds[0]);
	  int flds0len = strlen(flds[0]);
	  if (flds[0][flds0len-1] % 2 == 0)
	       ++flds[0][flds0len-1];
	  else
	       --flds[0][flds0len-1];
	  std::string readMate = std::string (flds[0]);
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

     if (stat (dirForGaps.c_str(), &statbuf) != 0) {
	  charb cmd(100);
	  sprintf (cmd, "mkdir %s", dirForGaps.c_str());
	  fprintf (stderr, "%s\n", (char *) cmd);
	  system (cmd);
     }
     int outputGroupNum = 0;
     while (1) {
	  bool isAtEnd = loadNeededReads();
	  if (readSeq.empty())
	       break;
	  charb outfileName(100);
	  sprintf (outfileName, "%s/readFile.%03d", dirForGaps.c_str(), outputGroupNum);
	  outfile.open (outfileName);
	  
	  for (unsigned int grp=0; grp<readsInGroup.size(); grp++) {
	       charb istr2(20), newDir(20);
	       outfile << '>' << grp << "\n";
	       if (outputDirProgress)
		    fprintf (stderr, "\rGrp %d", grp);
	       
	       std::map<std::string, int> willBeOutput;
	       willBeOutput.clear();
	       std::vector<std::string> readsToOutput;
	       readsToOutput.clear();
	       for (unsigned int readNum=0; readNum<readsInGroup[grp].size(); readNum++) {
		    std::string readName = readsInGroup[grp][readNum];
		    std::map<std::string, int>::iterator it = willBeOutput.find (readName);
		    if (it != willBeOutput.end())
			 continue;
		    readsToOutput.push_back (readName);
		    willBeOutput[readName] = 2; }
	       for (unsigned int readNum=0; readNum<mateReadsInGroup[grp].size(); readNum++) {
		    std::string readName = mateReadsInGroup[grp][readNum];
		    std::map<std::string, int>::iterator it = willBeOutput.find (readName);
		    if (it != willBeOutput.end())
			 continue;
		    willBeOutput[readName] = 1; }
	       std::set<std::string> alreadyOutput;
	       alreadyOutput.clear();
	       std::vector<std::string> bothMatesHaveMatches, mateBroughtInViaMatePair, readOnlys;
	       bothMatesHaveMatches.clear();
	       mateBroughtInViaMatePair.clear();
	       readOnlys.clear();
	       for (unsigned int readNum=0; readNum<readsToOutput.size(); readNum++) {
		    std::string readName = readsToOutput[readNum];
		    if (alreadyOutput.find(readName) != alreadyOutput.end())
			 continue;
		    std::string readMate = getReadMateName (readName);
		    alreadyOutput.insert (readName);
		    std::map<std::string, int>::iterator it2 = willBeOutput.find (readMate);
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
		    std::map<std::string,std::string>::iterator readSeqIt =
			 readSeq.find (bothMatesHaveMatches[readNum]);
		    if (readSeqIt == readSeq.end())
			 continue;
		    outfile << bothMatesHaveMatches[readNum] << " B\n" <<
			 readSeqIt->second << "\n";
	       }
	       int outCount = 0;
	       for (unsigned int readNum=0; readNum<mateBroughtInViaMatePair.size(); readNum++) {
		    ++outCount;
		    std::string extraStr;
		    std::string readName = mateBroughtInViaMatePair[readNum];
		    std::map<std::string,std::string>::iterator readSeqIt =
			 readSeq.find (readName);
		    if (readSeqIt == readSeq.end())
			 continue;
		    if (outCount % 2 == 1)
			 extraStr = std::string ("match");
		    else
			 extraStr = std::string ("mate");
		    outfile << readName << ' ' << extraStr << " M\n" <<
			 readSeq[readName] << "\n";
	       }
	       for (unsigned int readNum=0; readNum<readOnlys.size(); readNum++) {
		    std::string readName = readOnlys[readNum];
		    std::map<std::string,std::string>::iterator readSeqIt =
			 readSeq.find (readName);
		    if (readSeqIt == readSeq.end())
			 continue;
		    std::string mateRead = getReadMateName (readName);
		    outfile << readName << " O\n" << readSeq[readName] << '\n'
			    << mateRead << " N\n" << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";
	       }
	  }
	  outfile.close();
	  ++outputGroupNum;
	  if (outputDirProgress)
	       fprintf (stderr, "\n");
	       
	  if (isAtEnd)
	       break;
     }

     return (0);
}

bool loadNeededReads (void)
{
     FILE *infile;
     static unsigned int currentFileNum = 0;
     static off_t currentFileOffset = (off_t) 0;
     std::string localReadsFile;
     int64_t numReadSeqsLoaded;
     
     readSeq.clear();
     numReadSeqsLoaded = 0;
     while (currentFileNum < readFiles.size()) {
	  localReadsFile = readFiles[currentFileNum];
	  infile = fopen (localReadsFile.c_str(), "r");
	  fseeko (infile, currentFileOffset, SEEK_SET);
	  int offsetStart = localReadsFile.length()-5;
	  charb line(100), readNameStr(100);
	  char *cptr = line+1;
	  std::string readName;
	  bool isNeeded = false;
	  if (offsetStart < 0) offsetStart = 0;
	  if (localReadsFile.substr(offsetStart) == std::string ("fastq"))
	       goto fastQ;
	  while (fgets (line, 100, infile)) {
	       int lineLen;
	       cptr = line+1;
	       if (line[0] == '>') {
		    char *saveptr;
		    lineLen = strlen(line);
		    readNameStr = strtok_r(cptr, " \t\n", &saveptr);
		    readName = std::string (readNameStr);
		    if (readIsNeeded.find(readName) != readIsNeeded.end())
			 isNeeded = true;
		    else
			 isNeeded = false;
		    // What to do if we should stop before this read
		    if ((((numReadSeqsLoaded == maxReadsInMemory-1) && (readNameStr[strlen(readNameStr)-1] % 2 == 0)) ||
			 (numReadSeqsLoaded == maxReadsInMemory)) && isNeeded) {
			 currentFileOffset = ftello (infile);
			 currentFileOffset -= lineLen;
			 fclose (infile);
			 return false;
		    }
		    continue; }
	       if (isNeeded) {
		    line[strlen(line)-1] = 0; // chop (line);
		    if (readSeq.find(readName) == readSeq.end()) {
			 readSeq[readName] = std::string ("");
			 ++numReadSeqsLoaded; }
		    readSeq[readName] += std::string (line);
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
	       readName = std::string (readNameStr);
	       if (readIsNeeded.find(readName) != readIsNeeded.end())
		    isNeeded = true;
	       else
		    isNeeded = false;
	       // What to do if we should stop before this read
	       if ((((numReadSeqsLoaded == maxReadsInMemory-1) && (readNameStr[strlen(readNameStr)-1] % 2 == 0)) ||
		    (numReadSeqsLoaded == maxReadsInMemory)) && isNeeded) {
		    currentFileOffset = ftello (infile);
		    currentFileOffset -= lineLen;
		    fclose (infile);
		    return false;
	       }

	       fgets (line, 100, infile);
	       line[strlen(line)-1] = 0; // chop (line);
	       if (isNeeded) {
		    readSeq[readName] = std::string (line);
		    ++numReadSeqsLoaded;
	       }
	       fgets (line, 100, infile); // Skip quality bases
	       fgets (line, 100, infile);
	  }
//	  std::cout << "num read seqs loaded = " << numReadSeqsLoaded << std::endl;
     endProcessingReads:
	  fclose (infile);
	  ++currentFileNum;
	  currentFileOffset = (off_t) 0;
     }
     return true;
}

std::string getReadMateName (std::string readName)
{
     std::string readMate;
     int readNameLen = readName.length();
     charb readNameStr(100);
     if (readNameLen == 0)
	  return (std::string (""));
     strcpy (readNameStr, readName.c_str());
     if (readNameStr[readNameLen-1] % 2 == 0)
	       ++readNameStr[readNameLen-1];
	  else
	       --readNameStr[readNameLen-1];
     readMate = std::string (readNameStr);

     return (readMate);
}

void processArgs (int argc, char **argv)
{
     struct stat     statbuf;
     int fail = 0;
     outputDirProgress = false;

     maxReadsInMemory = 1000000000;
     for (int i=1; i<argc; i++) {
	  if (! strcmp(argv[i], "--faux-reads-file")) {
	       ++i;
	       fauxReadsFile = std::string(argv[i]);
	       continue; }
	  if (! strcmp(argv[i], "--faux-read-matches-to-kunis-file")) {
	       ++i;
	       fauxReadMatchesToKUnisFile = std::string(argv[i]);
	       continue; }
	  if (! strcmp(argv[i], "--read-matches-to-kunis-file")) {
	       ++i;
	       readMatchesToKUnisFile = std::string(argv[i]);
	       continue; }
	  if (! strcmp(argv[i], "--reads-file")) {
	       ++i;
	       readFiles.push_back (std::string (argv[i]));
	       continue; }
	  if (! strcmp(argv[i], "--dir-for-gaps")) {
	       ++i;
	       dirForGaps = std::string(argv[i]);
	       continue; }
	  if (! strcmp(argv[i], "--max-reads-to-allow")) {
	       ++i;
	       maxReadsInMemory = atoll (argv[i]);
	       continue; }
	  if (! strcmp(argv[i], "--output-dir-progress")) {
	       outputDirProgress = true;
	       continue; }
     }
     if (stat (fauxReadsFile.c_str(), &statbuf) != 0) {
	  std::cerr << "Faux reads file '" << fauxReadsFile << "' doesn't exist.\n";
	  fail = 1; }
     if (stat (fauxReadMatchesToKUnisFile.c_str(), &statbuf) != 0) {
	  std::cerr << "File of faux read matches to k-unitigs '" << fauxReadMatchesToKUnisFile << "' doesn't exist.\n";
	  fail = 1; }
     if (stat (readMatchesToKUnisFile.c_str(), &statbuf) != 0) {
	  std::cerr << "File of read matches to k-unitigs '" << readMatchesToKUnisFile << "' doesn't exist.\n";
	  fail = 1; }
     for (unsigned int i=0; i<readFiles.size(); i++) {
	  std::string readFile = readFiles[i];
	  if (stat (readFile.c_str(), &statbuf) != 0) {
	       std::cerr << "Read file '" << readFile << "' doesn't exist.\n";
	       fail = 1; }
     }
     
     if (fail)
	  exit (1);
}


     
