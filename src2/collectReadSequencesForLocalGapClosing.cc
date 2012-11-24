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
#include <signal.h>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <charbuf.hpp>
#include <thread_pool.hpp>
#include <pthread.h>
#include <err.hpp>
#include <misc.hpp>
#include <src2/collectReadSequencesForLocalGapClosing_cmdline.hpp>

struct arguments {
     int dirNum;
     charb fauxReadFileDataStr;
     basic_charb<remaper<char> > readFileDataStr;
     charb *superReadFastaString;
     charb *readPlacementString;
};

FILE *fopen_wait(int timeout_, char *filename){
int timeout=0;
FILE *infile;
do {
        infile = fopen (filename, "r");
        if(infile != NULL) break;
        sleep(1);
        timeout++;
        }
while(timeout < timeout_);
     if(infile == NULL)
          printf("Timed out on file %s\n", filename);
return(infile);
}

int reportNumGaps (const char *contigEndSeqFile);
int analyzeGap(struct arguments threadArg); // Returns int for now

cmdline_parse args;

std::string exeDir;

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
void loadNeededReads (setOfStrings &readIsNeeded, readNameToReadSequence &readSeq);
void checkArgs (void);


int main(int argc, char **argv)
{
     FILE *infile; // , *outfile;
     FILE *contigEndSeqFile;
     struct arguments threadArgs;
     charb line(100);
     vectorOfStrings fauxReadGroups;
     stringToIntMap group;
     std::unordered_map<stdString, char> end;
     setOfStrings readIsNeeded;
     readNameToReadSequence readSeq; // Read name to read sequence
     struct stat     statbuf;
     charb tempBuffer(100);
     charb *superReadFastaStrings, *readPlacementStrings;
     args.parse (argc, argv);
     checkArgs ();
     thread_pool<struct arguments, int> pool (args.num_threads_arg, analyzeGap);
     char *tempPtr = strrchr (argv[0], '/');
     if (tempPtr == NULL)
	  exeDir = std::string(".");
     else {
	  unsigned int diff = tempPtr - argv[0];
	  strcpy (tempBuffer, argv[0]);
	  tempBuffer[diff] = 0;
	  exeDir = std::string (tempBuffer); }

     contigEndSeqFile = fopen (args.contig_end_sequence_file_arg, "r");
     if (contigEndSeqFile == NULL) {
	  fprintf (stderr, "File '%s' doesn't exist! Bye!\n", args.contig_end_sequence_file_arg);
	  exit (1); }
     fclose (contigEndSeqFile);
     int numGaps = reportNumGaps (args.contig_end_sequence_file_arg);
     superReadFastaStrings = new charb[numGaps];
     readPlacementStrings = new charb[numGaps];

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
     
     contigEndSeqFile = fopen (args.contig_end_sequence_file_arg, "r");

     // Here we make the output directory
     if (stat (args.output_dir_arg, &statbuf) != 0) {
	  charb cmd(100);
	  sprintf (cmd, "mkdir %s", args.output_dir_arg);
	  system (cmd); }

//     int outputGroupNum = 0;
#define AFTER
#ifdef AFTER
     charb tempLine(1000);

     loadNeededReads (readIsNeeded, readSeq);
     for (unsigned int dirNum=0; dirNum<readsInGroup.size(); dirNum++) {
	  charb istr2(20), newDir(20);
	  threadArgs.fauxReadFileDataStr.clear();
	  threadArgs.readFileDataStr.clear();
	  for (int j=0; j<4; j++) {
	       fgets (tempLine, 1000, contigEndSeqFile);
	       strcat (threadArgs.fauxReadFileDataStr, tempLine); }
	  stringToIntMap willBeOutput;
	  willBeOutput.clear();
	  vectorOfStrings readsToOutput;
	  readsToOutput.clear();
	  for (unsigned int readNum=0; readNum<readsInGroup[dirNum].size(); readNum++) {
	       stdString readName = readsInGroup[dirNum][readNum];
	       stringToIntMap::iterator it = willBeOutput.find (readName);
	       if (it != willBeOutput.end())
		    continue;
	       readsToOutput.push_back (readName);
	       willBeOutput[readName] = 2; }
	  for (unsigned int readNum=0; readNum<mateReadsInGroup[dirNum].size(); readNum++) {
	       stdString readName = mateReadsInGroup[dirNum][readNum];
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
	       sprintf (line, ">%s B\n%s\n", bothMatesHaveMatches[readNum].c_str(), (readSeqIt->second).c_str());
	       strcat (threadArgs.readFileDataStr, line);
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
	       sprintf (line, ">%s %s M\n%s\n", readName.c_str(), extraStr.c_str(), (readSeq[readName]).c_str());
	       strcat (threadArgs.readFileDataStr, line);
	  }

	  int tempInt;
	  for (unsigned int readNum=0; readNum<readOnlys.size(); readNum++) {
	       stdString readName = readOnlys[readNum];
	       stringToStringMap::iterator readSeqIt =
		    readSeq.find (readName);
	       if (readSeqIt == readSeq.end())
		    continue;
	       stdString mateRead = getReadMateName (readName);
	       sprintf (line, ">%s O\n%s\n>%s N\nN\n", readName.c_str(), (readSeq[readName]).c_str(), mateRead.c_str());
	       strcat (threadArgs.readFileDataStr, line);
	  }
	  threadArgs.dirNum = dirNum;
	  threadArgs.superReadFastaString = &superReadFastaStrings[dirNum];
          threadArgs.readPlacementString = &readPlacementStrings[dirNum];
	  // Get next available thread when available
	  pool.submit_job (&threadArgs, &tempInt);
	  
     }
     
     pool.release_workers();

     charb outfileName;
     sprintf (outfileName, "./superReadSequences.fasta"); // , args.output_dir_arg);
     FILE *outfile = fopen ((char *)outfileName, "w");
     for (int i=0; i<numGaps; i++) {
	  if (superReadFastaStrings[i] == NULL)
	       continue;
	  fputs (superReadFastaStrings[i], outfile); }
     fclose (outfile);
     
     sprintf (outfileName, "./readPlacementsInSuperReads.final.read.superRead.offset.ori.txt"); // , args.output_dir_arg);
     outfile = fopen ((char *)outfileName, "w");
     for (int i=0; i<numGaps; i++) {
	  if (readPlacementStrings[i] == NULL)
	       continue;
	  fputs (readPlacementStrings[i], outfile); }
     fclose (outfile);
     
     return 0;
}




#else
     while (1) {
//	  bool isAtEnd = loadNeededReads();
	  loadNeededReads (readIsNeeded, readSeq);
	  if (readSeq.empty())
	       break;
	  charb outfileName(100);
//	  sprintf (outfileName, "%s/readFile.%03d", args.dir_for_gaps_arg, outputGroupNum);
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
//	  ++outputGroupNum;
	       
	  if (isAtEnd)
	       break;
     }

     return (0);
}
#endif

void loadNeededReads (setOfStrings &readIsNeeded, readNameToReadSequence &readSeq)
{
     FILE *infile;
     static unsigned int currentFileNum = 0;
//     static off_t currentFileOffset = (off_t) 0;
     const char *localReadsFile;
     uint64_t numReadSeqsLoaded;
     
     readSeq.clear();
     numReadSeqsLoaded = 0;
     while (currentFileNum < args.reads_file_arg.size()) {
	  localReadsFile = args.reads_file_arg[currentFileNum];
	  infile = fopen (localReadsFile, "r");
//	  fseeko (infile, currentFileOffset, SEEK_SET);
	  int offsetStart = strlen (localReadsFile)-5;
	  charb line(100), readNameStr(100);
	  char *cptr = line+1;
	  stdString readName;
	  bool isNeeded = false;
	  if (offsetStart < 0) offsetStart = 0;
	  const char *cptr2 = localReadsFile + offsetStart;
	  if (strcmp (cptr2, "fastq") == 0)
	       goto fastQ;
	  // If we get here it's a fasta file
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
//	  currentFileOffset = (off_t) 0;
     }
     return;
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

int analyzeGap(struct arguments threadArg)
{
     struct stat statbuf;
     FILE *infile, * outfile;

     // Doing the actual work of the worker thread
     charb outDirName(100), tempFileName(10), cmd(100), line(100);

     sprintf (outDirName, "%s/gap%09ddir", args.output_dir_arg, threadArg.dirNum);
    
     if (stat (outDirName, &statbuf) != 0) {
	  sprintf (cmd, "mkdir %s", (char *)outDirName);
	  system (cmd); }
     sprintf (tempFileName, "%s/fauxReads.fasta", (char *) outDirName);
     outfile = fopen (tempFileName, "w");
     fputs (threadArg.fauxReadFileDataStr, outfile);
     fclose (outfile);
     sprintf (tempFileName, "%s/reads.fasta", (char *) outDirName);
     outfile = fopen (tempFileName, "w");
     fputs (threadArg.readFileDataStr, outfile);
     fclose (outfile);
     sprintf (cmd, "%s/closeGaps.oneDirectory.perl --dir-to-change-to %s --Celera-terminator-directory %s --reads-file reads.fasta --output-directory outputDir --max-kmer-len %d --min-kmer-len %d --maxnodes %d --mean-for-faux-inserts %d --stdev-for-faux-inserts %d --use-all-kunitigs --noclean 1>%s/out.err 2>&1", exeDir.c_str(), (char *) outDirName, args.Celera_terminator_directory_arg, args.max_kmer_len_arg, args.min_kmer_len_arg, args.max_nodes_arg, args.mean_for_faux_inserts_arg, args.stdev_for_faux_inserts_arg, (char *) outDirName, (char *) outDirName);
     printf ("Working on dir %s on thread %ld\n", (char *) outDirName, pthread_self());
     sprintf (tempFileName, "%s/passingKMer.txt", (char *) outDirName);
     int passingKMerValue = 0;
     system (cmd);
     /* Now, if "passingKMer.txt" exists in outDirName, copy the files
	superReadSequences.fasta and
	readPlacementsInSuperReads.final.read.superRead.offset.ori.txt (after appropriate
	modifications), copy to the desired output directory from (e.g.)
	'outDirName'/work_localReadsFile_41_2 */

     while(passingKMerValue <11){
     	infile = fopen_wait (60,tempFileName);
	fgets(line,100,infile);
	passingKMerValue=atoi(line);
     	fclose (infile);
     }

     if(passingKMerValue == 11){
          if (! args.keep_directories_flag){
               sprintf (cmd, "rm -rf %s", (char *) outDirName);
               system ((char *) cmd);
		}
          return (0);
        }
     
     // If we get here we have found a join and passingKMer.txt exists
     sprintf (tempFileName, "%s/work_localReadsFile_%d_2/superReadSequences.fasta", (char *) outDirName, passingKMerValue);
     infile = fopen_wait (10,tempFileName);
     if(infile == NULL){
          if (! args.keep_directories_flag)
               sprintf (cmd, "rm -rf %s", (char *) outDirName);
          system ((char *) cmd);
          return (0);
        }        
     fgets (line, 100, infile);
     sprintf (*(threadArg.superReadFastaString), ">%d\n", threadArg.dirNum);
     while (fgets (line, 100, infile))
	  strcat (*(threadArg.superReadFastaString), line);
//     while (fgets_append (*(threadArg.superReadFastaString), infile))
//          ;
     fclose (infile);
     sprintf (tempFileName, "%s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt", (char *) outDirName, passingKMerValue);
        infile = fopen_wait (10,tempFileName);
     if(infile == NULL){
          if (! args.keep_directories_flag)
               sprintf (cmd, "rm -rf %s", (char *) outDirName);
          system ((char *) cmd);
          return (0);
        }
     ExpBuffer<char *>flds;
     charb readPlacementLine;
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
          sprintf (readPlacementLine, "%s %d %s %s\n", flds[0], threadArg.dirNum, flds[2], flds[3]);
          strcat (*(threadArg.readPlacementString), readPlacementLine);
     }
     fclose (infile);
#if 0
     sprintf (cmd, "cp %s/work_localReadsFile_%d_2/superReadSequences.fasta %s/superReadSequences.%09d.fasta", (char *) outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
     system (cmd);
     sprintf (cmd, "cp %s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt %s/readPlacementsInSuperReads.final.read.superRead.offset.ori.%09d.txt", (char *) outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
     system (cmd);
#endif

//printf("Gap %s passing kmer %d super read %s\n",(char *) outDirName,passingKMerValue,(char*)*(threadArg.superReadFastaString));
     if (! args.keep_directories_flag) {
	  sprintf (cmd, "rm -rf %s", (char *) outDirName);
	  system ((char *) cmd);
     }

     return (0);
}

int reportNumGaps (const char *fn)
{
     charb cmd(100), line(100);
     FILE *infile;

     sprintf (cmd, "tail -2 %s | head -1", fn);
     infile = popen (cmd, "r");
     fgets (line, 100, infile);
     char *cptr = line;
     while (! isdigit(*cptr))
          cptr++;
     int lastFauxContig = atoi (cptr);
     int numGaps = (lastFauxContig+1)/2;
     return (numGaps);
}

