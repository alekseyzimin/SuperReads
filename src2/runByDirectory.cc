#include <sys/types.h>
#include <sys/wait.h>
#include <ftw.h>
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
       
int remove_function(const char *fpath, const struct stat *sb,
                          int typeflag, struct FTW *ftwbuf){
return(remove(fpath));
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
     args.parse (argc, argv);
     //thread_pool<struct arguments, int> pool (args.num_threads_arg, analyzeGap);

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
     if (stat (args.output_dir_arg, &statbuf) != 0)
          mkdir((char *)args.output_dir_arg, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

     pid_t pid;
     int status=0;
     int dirNum;

     for (dirNum=0; ! atEnd; dirNum++) {
     if(dirNum>=args.num_threads_arg)
             wait(&status);

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
	  // Get next available thread when available
	  //pool.submit_job (&threadArgs, &tempInt);
          fflush(stdout);
          pid=fork();
          if(pid<0)
                perror("Fork failed");
          else if(pid==0){
                        analyzeGap(threadArgs);
                        exit(0);
                }
     }

        for(int ttt=0;ttt<args.num_threads_arg;ttt++)
        	wait(&status);
        

     //pool.release_workers();

        if (! args.keep_directories_flag)
            nftw((char *) args.output_dir_arg,remove_function,2048,FTW_DEPTH|FTW_PHYS);

     return 0;
}

int analyzeGap(struct arguments threadArg)
{
     struct stat statbuf;
     FILE *infile, * outfile;

     // Doing the actual work of the worker thread
     charb outDirName(100), tempFileName(10), cmd(100), line(100);

     sprintf (outDirName, "%s/gap%09ddir", args.output_dir_arg, threadArg.dirNum);
     if (stat ((char *)outDirName, &statbuf) != 0)
          mkdir((char *)outDirName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 
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
        if (infile == NULL) {
             passingKMerValue = 11;
             break; }
        fgets(line,100,infile);
        passingKMerValue=atoi(line);
        fclose (infile);
     }

     if(passingKMerValue == 11){
          if (! args.keep_directories_flag)
               nftw((char *) outDirName,remove_function,2048,FTW_DEPTH|FTW_PHYS);
          return (0);
        }
     charb superReadFastaString(1000);
     sprintf (tempFileName, "%s/work_localReadsFile_%d_2/superReadSequences.fasta", (char *) outDirName, passingKMerValue);
     infile = fopen_wait (10,tempFileName);
     if(infile == NULL){
          fprintf(stderr,"failed to open superReadSequences.fasta %s passingKMerValue %d\n", (char *) outDirName,passingKMerValue);
          if (! args.keep_directories_flag)
               nftw((char *) outDirName,remove_function,2048,FTW_DEPTH|FTW_PHYS);
          return (0);
        }

     fgets (line, 100, infile);
     sprintf (superReadFastaString, "%d ", threadArg.dirNum);
     while (fgets (line, 100, infile)){
          line.chomp();
          strcat (superReadFastaString, line);
        }
     fclose (infile);
     sprintf (tempFileName, "%s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt", (char *) outDirName, passingKMerValue);
        infile = fopen_wait (10,tempFileName);
     if(infile == NULL){
          fprintf(stderr,"failed to open readPlacementsInSuperReads.final.read.superRead.offset.ori.txt %s passingKMerValue %d\n", (char *) outDirName,passingKMerValue);
          if (! args.keep_directories_flag)
               nftw((char *) outDirName,remove_function,2048,FTW_DEPTH|FTW_PHYS);
          return (0);
        }
     ExpBuffer<char *>flds;
     charb readPlacementLine,readPlacementLines;
     while (fgets (line, 100, infile)) {
          getFldsFromLine (line, flds);
          sprintf (readPlacementLine, " %s %d %s %s", flds[0], threadArg.dirNum, flds[2], flds[3]);
          strcat(readPlacementLines,readPlacementLine);
     }

     fclose (infile);
     fprintf(stderr,"%s%s\n",(char*)superReadFastaString,(char*)readPlacementLines);

#if 0
     sprintf (cmd, "cp %s/work_localReadsFile_%d_2/superReadSequences.fasta %s/superReadSequences.%09d.fasta", (char *) outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
     system (cmd);
     sprintf (cmd, "cp %s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt %s/readPlacementsInSuperReads.final.read.superRead.offset.ori.%09d.txt", (char *) outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
     system (cmd);
#endif

//printf("Gap %s passing kmer %d super read %s\n",(char *) outDirName,passingKMerValue,(char*)*(threadArg.superReadFastaString));
     if (! args.keep_directories_flag)
            nftw((char *) outDirName,remove_function,2048,FTW_DEPTH|FTW_PHYS);

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

