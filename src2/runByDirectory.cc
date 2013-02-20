#include <sys/types.h>
#include <sys/wait.h>
#include <ftw.h>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <string>
#include <sys/stat.h>
#include <pthread.h>

#include <vector>
#include <iostream>
#include <sstream>

#include <charbuf.hpp>
#include <thread_pool.hpp>
#include <misc.hpp>
#include <jellyfish/err.hpp>
#include <src2/runByDirectory_cmdline.hpp>

struct arguments {
  float mean;
  float stdev;
  int   dirNum;
  charb fauxReadFileDataStr;
  basic_charb<remaper<char> > readFileDataStr;

  void clear() {
    fauxReadFileDataStr.clear();
    readFileDataStr.clear();
  }
};

// Version of standard calls which throw a runtime_error in case of an
// error.
FILE* fopen_throw(const char* path, const char* mode, const char* msg = "") {
  FILE* res = fopen(path, mode);
  if(res == 0)
    eraise(std::runtime_error) << msg
                               << "\nfopen(" << path << ", " << mode << ") failed: "
                               << jellyfish::err::no;
  return res;
}

void mkdir_throw(const char* path, mode_t mode, const char* msg = "") {
  int res = mkdir(path, mode);
  if(res == -1)
    eraise(std::runtime_error) << msg
                               << "\nmkdir(" << path << ") failed: "
                               << jellyfish::err::no;
}

void system_throw(const char* command, const char* msg = "") {
  int res = system(command);

  switch(res) {
  case 0: return;
  case -1: eraise(std::runtime_error) << msg
                                      << "\nsystem(" << command << ") failed: "
                                      << jellyfish::err::no;
  default:
    break;
  }
  // If get there, system call succeeded but the command itself
  // failed.
  if(WIFEXITED(res))
    eraise(std::runtime_error) << msg
                               << "\ncommand '" << command << "' returned error code "
                               << WEXITSTATUS(res);
  if(WIFSIGNALED(res))
    eraise(std::runtime_error) << msg
                               << "\ncommand '" << command << "' killed by signal "
                               << strsignal(WTERMSIG(res));
  eraise(std::runtime_error) << msg
                             << "\ncommand '" << command 
                             << "' terminated for an unknown reason";
}

// FILE *fopen_wait(int timeout_, char *filename) {
//   int   timeout  = 0;
//   FILE *infile;
//   do {
//     infile = fopen (filename, "r");
//     if(infile != NULL) break;
//     sleep(1);
//     timeout++;
//   } while(timeout < timeout_);
//   if(infile == NULL)
//     printf("Timed out on file %s\n", filename);
//   return(infile);
// }

int remove_function(const char *fpath, const struct stat *sb,
                    int typeflag, struct FTW *ftwbuf) {
  return remove(fpath);
}
/// Equivalent to 'rm -rf path'
void rm_rf(const char* path) {
  nftw(path, remove_function, 2048, FTW_DEPTH|FTW_PHYS);
}

/// Attempt to set the close on exec flag for the file descriptor. Any
/// error is ignored.
void close_on_exec(FILE* file) {
  if(!file)
    return;
  int fd = fileno(file);
  if(fd == -1)
    die << "fileno" << jellyfish::err::no;
  int flags = fcntl(fd, F_GETFD);
  if(flags == -1)
    die << "fcntl F_GETFD" << jellyfish::err::no;
  flags = fcntl(fd, F_SETFD, flags|FD_CLOEXEC);
  if(flags == -1)
    die << "fcntl F_SETFD" << jellyfish::err::no;
}

int reportNumGaps (const char *contigEndSeqFile);
int analyzeGap(struct arguments threadArg); // Returns int for now

// THE NEXT 2 LINES ARE COPIED FROM GUILLAUME
// GLOBAL: command line switches
cmdline_parse args;

std::string exeDir;

int main (int argc, char **argv)
{
  args.parse (argc, argv);

  struct arguments    threadArgs;
  FILE               *contigEndSeqFile;
  FILE               *meanAndStdevFile;
  std::vector<FILE*>  readSeqFilesByDir;
  charb               tempBuffer(100);

  char *tempPtr = strrchr (argv[0], '/');
  if(tempPtr == NULL) {
    exeDir = std::string(".");
  } else {
    exeDir = std::string(argv[0], tempPtr - argv[0]);
  }

  contigEndSeqFile = fopen(args.contig_end_sequence_file_arg, "r");
  close_on_exec(contigEndSeqFile);
  if(!contigEndSeqFile)
    die << "Can't open '" << args.contig_end_sequence_file_arg << "'"
        << jellyfish::err::no;

  charb fn(200), line(100), tempLine(1000);

  for (int fileNum=0; 1; fileNum++) {
    sprintf (fn, "%s/readFile.%03d", args.dir_for_read_sequences_arg, fileNum);
    FILE *filetemp = fopen(fn, "r");
    if(!filetemp)
      continue;
    close_on_exec(filetemp);
    readSeqFilesByDir.push_back (filetemp);
    fgets(line, 100, filetemp); // Gets past the first line ('>0')
  }


  // Here we make the output directory
  int res = mkdir(args.output_dir_arg, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(res == -1 && errno != EEXIST)
    die << "Can't create output directory '" << args.output_dir_arg
        << jellyfish::err::no;


  meanAndStdevFile = NULL;
  if(args.mean_and_stdev_file_given) {
    meanAndStdevFile = fopen(args.mean_and_stdev_file_arg, "r");
    close_on_exec(meanAndStdevFile);
    if(!meanAndStdevFile)
      die << "Can't open file '" << args.mean_and_stdev_file_arg << "'"
          << jellyfish::err::no;
  }

  int   status = 0;
  int   dirNum;
  bool  atEnd  = false;
  for (dirNum=0; ! atEnd; dirNum++) {
    if(dirNum >= args.num_threads_arg)
      wait(&status);

    threadArgs.clear();
    for (int j = 0; j < 4; ++j) {
      fgets (tempLine, 1000, contigEndSeqFile);
      strcat (threadArgs.fauxReadFileDataStr, tempLine);
    }
    //	       fgets_append (threadArgs.fauxReadFileDataStr, contigEndSeqFile); }

    atEnd = true;
    for (unsigned int i = 0; i < readSeqFilesByDir.size(); ++i) {
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

    if (args.mean_and_stdev_file_given) {
      int dirNumTemp;
      int scanned = fscanf (meanAndStdevFile, "%d %f %f\n",
                            &dirNumTemp, &threadArgs.mean, &threadArgs.stdev);
      if(scanned != 3)
        die << "Failed to parse file '" << args.mean_and_stdev_file_arg << "'"
            << jellyfish::err::no;
    } else {
      threadArgs.mean  = args.mean_for_faux_inserts_arg;
      threadArgs.stdev = args.stdev_for_faux_inserts_arg;
    }
    threadArgs.dirNum = dirNum;

    // Get next available thread when available
    //pool.submit_job (&threadArgs, &tempInt);
    fflush(stdout);
    switch(fork()) {
    case -1:
      die << "Fork failed" << jellyfish::err::no;
    case 0:
      analyzeGap(threadArgs);
      exit(0);
    default:
      break;
    }
  }

  for(int ttt=0;ttt<args.num_threads_arg;ttt++)
    wait(&status);

  if(!args.keep_directories_flag)
    rm_rf(args.output_dir_arg);

  return 0;
}

// Doing the actual work of the worker thread
void do_analyzeGap(struct arguments& threadArg, const char* outDirName);
int analyzeGap(struct arguments threadArg) {
  charb outDirName(100);

  sprintf (outDirName, "%s/gap%09ddir", args.output_dir_arg, threadArg.dirNum);
  int res = mkdir((char *)outDirName, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
  if(res == -1 && errno != EEXIST) // OK if directory already exists
    return -1;

  try {
    do_analyzeGap(threadArg, outDirName);
  } catch(std::runtime_error e) {
    fprintf(stderr, "Analyze gap failed for dir '%s': %s\n",
            (const char*)outDirName, e.what());
    if(!args.keep_directories_flag)
      rm_rf(outDirName);
    return -1;
  }
  return 0;
}

void do_analyzeGap(struct arguments& threadArg, const char* outDirName) {
  FILE *infile, * outfile;

  charb tempFileName(10), line(100);

  sprintf (tempFileName, "%s/fauxReads.fasta", outDirName);
  outfile = fopen_throw (tempFileName, "w", "Write faux read sequence");
  fputs (threadArg.fauxReadFileDataStr, outfile);
  fclose (outfile);

  sprintf (tempFileName, "%s/reads.fasta", outDirName);
  outfile = fopen_throw (tempFileName, "w", "Write read sequence");
  fputs (threadArg.readFileDataStr, outfile);
  fclose (outfile);

  std::ostringstream cmd;
  cmd << exeDir.c_str() << "/closeGaps.oneDirectory.perl"
      << " 1>" << outDirName << "/out.err 2>&1"
      << " --dir-to-change-to " << outDirName
      << " --Celera-terminator-directory " << args.Celera_terminator_directory_arg
      << " --reads-file reads.fasta --output-directory outputDir"
      << " --max-kmer-len " << args.max_kmer_len_arg
      << " --min-kmer-len " << args.min_kmer_len_arg
      << " --maxnodes " << args.max_nodes_arg
      << " --mean-for-faux-inserts " << threadArg.mean
      << " --stdev-for-faux-inserts " << threadArg.stdev
      << " --num-stdevs-allowed " << args.num_stdevs_allowed_arg
      << " --use-all-kunitigs --noclean";

  printf ("Working on dir %s on thread %ld\n", outDirName, pthread_self());
  sprintf (tempFileName, "%s/passingKMer.txt", outDirName);
  system_throw(cmd.str().c_str());

  /* Now, if "passingKMer.txt" exists in outDirName, copy the files
     superReadSequences.fasta and
     readPlacementsInSuperReads.final.read.superRead.offset.ori.txt (after appropriate
     modifications), copy to the desired output directory from (e.g.)
     'outDirName'/work_localReadsFile_41_2 */
  int passingKMerValue = 0;

  // while(passingKMerValue <11){
  //   infile = fopen_wait (60,tempFileName);
  //   if (infile == NULL) {
  //     passingKMerValue = 11;
  //     break; }
  //   fgets(line,100,infile);
  //   passingKMerValue=atoi(line);
  //   fclose (infile);
  // }

  // if(passingKMerValue == 11)
  //   return;

  infile = fopen_throw(tempFileName, "r", "Reading passing k-mer size");
  int scanned = fscanf(infile, "%d", &passingKMerValue);
  if(scanned != 1)
    eraise(std::runtime_error) << "Failed to read passing k-mer size"
                               << jellyfish::err::no;

  charb superReadFastaString(1000);
  sprintf (tempFileName, "%s/work_localReadsFile_%d_2/superReadSequences.fasta",
           outDirName, passingKMerValue);
  infile = fopen_throw (tempFileName, "r", "Reading SuperRead sequences");
  fgets (line, 100, infile);
  sprintf (superReadFastaString, "%d ", threadArg.dirNum);
  while (fgets (line, 100, infile)){
    line.chomp();
    strcat (superReadFastaString, line);
  }
  fclose (infile);

  sprintf (tempFileName, "%s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt",
           outDirName, passingKMerValue);
  infile = fopen_throw (tempFileName, "r", "Reading read placements");
  ExpBuffer<char*> flds;
  charb readPlacementLine,readPlacementLines;
  while (fgets (line, 100, infile)) {
    getFldsFromLine (line, flds);
    sprintf(readPlacementLine, " %s %d %s %s", flds[0], threadArg.dirNum, flds[2], flds[3]);
    strcat(readPlacementLines,readPlacementLine);
  }
  fclose (infile);

  fprintf(stderr,"%s%s\n",(char*)superReadFastaString,(char*)readPlacementLines);

// #if 0
//   sprintf (cmd, "cp %s/work_localReadsFile_%d_2/superReadSequences.fasta %s/superReadSequences.%09d.fasta", outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
//   system (cmd);
//   sprintf (cmd, "cp %s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt %s/readPlacementsInSuperReads.final.read.superRead.offset.ori.%09d.txt", outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
//   system (cmd);
// #endif

  return;
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

