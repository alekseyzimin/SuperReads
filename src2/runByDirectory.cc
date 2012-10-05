#include <cstdlib>
#include <cstdio>
#include <signal.h>
#include <string>
#include <sys/stat.h>
#include <charbuf.hpp>
#include <thread_pool.hpp>
#include <pthread.h>
#include <vector>

struct arguments {
     int dirNum;
     charb fauxReadFileDataStr;
     charb readFileDataStr;
};

int analyzeGap(struct arguments threadArg); // Returns int for now
std::string outputDir = "subdir2";
std::string exeDir;
std::string CeleraTerminatorDirectory;
int maxKMerLen, minKMerLen;

#define NUM_THREADS 4

int main (int argc, char **argv)
{
     struct arguments threadArgs;
     std::string dirForContigEnds = ".";
     std::string dirForReadSeqs = "subdir";
     FILE *contigEndSeqFile;
     std::vector <FILE *> readSeqFilesByDir;
     struct stat statbuf;
     charb tempBuffer(100);
     thread_pool<struct arguments, int> pool (NUM_THREADS, analyzeGap);
     CeleraTerminatorDirectory = std::string ("/genome2/raid/tri/localGapClosing/rhodobacter/subdir");
     maxKMerLen = 65; minKMerLen = 17;

     char *tempPtr = strrchr (argv[0], '/');
     if (tempPtr == NULL)
	  exeDir = std::string(".");
     else {
	  unsigned int diff = tempPtr - argv[0];
	  strcpy (tempBuffer, argv[0]);
	  tempBuffer[diff] = 0;
	  exeDir = std::string (tempBuffer); }

     std::string contigEndSeqFilename = dirForContigEnds + "/contig_end_pairs.100.fa";
     contigEndSeqFile = fopen (contigEndSeqFilename.c_str(), "r");
     charb fn(200), line(100);

     for (int fileNum=0; 1; fileNum++) {
	  sprintf (fn, "%s/readFile.%03d", dirForReadSeqs.c_str(), fileNum);
	  if (stat (fn, &statbuf) != 0)
	       break;
	  FILE *filetemp = fopen (fn, "r");
	  readSeqFilesByDir.push_back (filetemp);
	  fgets (line, 100, filetemp); // Gets past the first line ('>0')
     }
     bool atEnd = false;

     // Here we make the output directory
     if (stat (outputDir.c_str(), &statbuf) != 0) {
	  charb cmd(100);
	  sprintf (cmd, "mkdir %s", outputDir.c_str());
	  system (cmd); }

     int dirNum;
     for (dirNum=0; ! atEnd; dirNum++) {
	  threadArgs.fauxReadFileDataStr.clear();
	  for (int j=0; j<4; j++) {
	       fgets_append (threadArgs.fauxReadFileDataStr, contigEndSeqFile); }
	  atEnd = true;
	  threadArgs.readFileDataStr.clear();
	  for (unsigned int i=0; i<readSeqFilesByDir.size(); i++) {
	       int lineNum=0;
	       while (fgets (line, 100, readSeqFilesByDir[i])) {
		    if (line[0] == '>') {
			 atEnd = false;
			 break;
		    }
		    if (lineNum % 2 == 0)
			 strcat (threadArgs.readFileDataStr, ">");
		    strcat (threadArgs.readFileDataStr, line);
		    ++lineNum;
	       }
	  }
	  threadArgs.dirNum=dirNum;
	  // Get next available thread when available
	  int tempInt;
	  pool.submit_job (&threadArgs, &tempInt);

     }
     
     pool.release_workers();
     return 0;
}

int analyzeGap(struct arguments threadArg)
{
     struct stat statbuf;
     FILE *infile;

     // Doing the actual work of the worker thread
     charb outDirName(100), outContigEndSeqFileName(100), outReadsInBucketFileName(100), cmd(100);

     sprintf (outDirName, "%s/gap%09ddir", outputDir.c_str(), threadArg.dirNum);
     if (stat (outDirName, &statbuf) != 0) {
	  sprintf (cmd, "mkdir %s", (char *)outDirName);
	  system (cmd); }
     sprintf (outContigEndSeqFileName, "%s/fauxReads.fasta", (char *) outDirName);
     FILE *outfile = fopen (outContigEndSeqFileName, "w");
     fputs (threadArg.fauxReadFileDataStr, outfile);
     fclose (outfile);
     sprintf (outReadsInBucketFileName, "%s/reads.fasta", (char *) outDirName);
     outfile = fopen (outReadsInBucketFileName, "w");
     fputs (threadArg.readFileDataStr, outfile);
     fclose (outfile);
     sprintf (cmd, "%s/closeGaps.oneDirectory.perl --dir-to-change-to %s --Celera-terminator-directory %s --reads-file reads.fasta --output-directory outputDir --max-kmer-len %d --min-kmer-len %d --use-all-kunitigs --noclean >%s/out.err 2>%s/out.err", exeDir.c_str(), (char *) outDirName, CeleraTerminatorDirectory.c_str(), maxKMerLen, minKMerLen, (char *) outDirName, (char *) outDirName);
     printf ("Working on dir %s\n", (char *) outDirName);
     system (cmd);
     /* Now, if "passingKMer.txt" exists in outDirName, copy the files
	superReadSequences.fasta and
	readPlacementsInSuperReads.final.read.superRead.offset.ori.txt (after appropriate
	modifications), copy to the desired output directory from (e.g.)
	'outDirName'/work_localReadsFile_41_2 */
     charb passingKMerFilename(100);
     sprintf (passingKMerFilename, "%s/passingKMer.txt", (char *) outDirName);
     if (stat (passingKMerFilename, &statbuf) == 0) { // It exists
	  charb goodWorkDirectory(100);
	  int passingKMerValue;
	  infile = fopen (passingKMerFilename, "r");
	  fscanf (infile, "%d", &passingKMerValue);
	  fclose (infile);
	  sprintf (cmd, "cp %s/work_localReadsFile_%d_2/superReadSequences.fasta %s/superReadSequences.%09d.fasta", (char *) outDirName, passingKMerValue, outputDir.c_str(), threadArg.dirNum);
	  system (cmd);
	  sprintf (cmd, "cp %s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt %s/readPlacementsInSuperReads.final.read.superRead.offset.ori.%09d.txt", (char *) outDirName, passingKMerValue, outputDir.c_str(), threadArg.dirNum);
	  system (cmd);
     }

     return (0);
}
