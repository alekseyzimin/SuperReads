#include <cstdlib>
#include <cstdio>
#include <signal.h>
#include <string>
#include <sys/stat.h>
#include <charbuf.hpp>
#include <thread_pool.hpp>
#include <misc.hpp>
#include <pthread.h>
#include <vector>
#include <src2/runByDirectory_cmdline.hpp>

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

// THE NEXT 2 LINES ARE COPIED FROM GUILLAUME
// GLOBAL: command line switches
cmdline_parse args;

std::string exeDir;

int main (int argc, char **argv)
{
     struct arguments threadArgs;
     FILE *contigEndSeqFile;
     std::vector <FILE *> readSeqFilesByDir;
     struct stat statbuf;
     charb tempBuffer(100);
     charb *superReadFastaStrings, *readPlacementStrings;
     args.parse (argc, argv);
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

     contigEndSeqFile = fopen (args.contig_end_sequence_file_arg, "r");
     charb fn(200), line(100), tempLine(1000);

     for (int fileNum=0; 1; fileNum++) {
	  sprintf (fn, "%s/readFile.%03d", args.dir_for_read_sequences_arg, fileNum);
	  if (stat (fn, &statbuf) != 0)
	       break;
	  FILE *filetemp = fopen (fn, "r");
	  readSeqFilesByDir.push_back (filetemp);
	  fgets (line, 100, filetemp); // Gets past the first line ('>0')
     }
     bool atEnd = false;

     // Here we make the output directory
     if (stat (args.output_dir_arg, &statbuf) != 0) {
	  charb cmd(100);
	  sprintf (cmd, "mkdir %s", args.output_dir_arg);
	  system (cmd); }

     int dirNum;
     int tempInt;
     for (dirNum=0; ! atEnd; dirNum++) {
	  threadArgs.fauxReadFileDataStr.clear();
	  for (int j=0; j<4; j++) {
	       fgets (tempLine, 1000, contigEndSeqFile);
	       strcat (threadArgs.fauxReadFileDataStr, tempLine); }
//	       fgets_append (threadArgs.fauxReadFileDataStr, contigEndSeqFile); }
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

