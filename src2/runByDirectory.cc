#include <sys/types.h>
#include <sys/wait.h>
#include <sys/file.h>
#include <ftw.h>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <string>
#include <sys/stat.h>

#include <vector>
#include <iostream>
#include <sstream>

#include <charbuf.hpp>
#include <thread_pool.hpp>
#include <misc.hpp>
#include <jellyfish/err.hpp>
#include <src2/runByDirectory_cmdline.hpp>

namespace err = jellyfish::err;

struct arguments {
  int mean;
  int stdev;
  long int   dirNum;
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
    throw std::runtime_error(err::msg() << msg
                             << "\nfopen(" << path << ", " << mode << ") failed: "
                             << err::no);
  return res;
}

void mkdir_throw(const char* path, mode_t mode, const char* msg = "") {
  int res = mkdir(path, mode);
  if(res == -1)
    throw std::runtime_error(err::msg() << msg
                             << "\nmkdir(" << path << ") failed: "
                             << ::err::no);
}

void system_throw(const char* command, const char* msg = "") {
  int res = system(command);

  switch(res) {
  case 0: return;
  case -1: throw std::runtime_error(err::msg() << msg
                                    << "\nsystem(" << command << ") failed: "
                                    << err::no);
  default:
    break;
  }
  // If get there, system call succeeded but the command itself
  // failed.
  if(WIFEXITED(res))
    throw std::runtime_error(err::msg() << msg
                             << "\ncommand '" << command << "' returned error code "
                             << WEXITSTATUS(res));
  if(WIFSIGNALED(res))
    throw std::runtime_error(err::msg() << msg
                             << "\ncommand '" << command << "' killed by signal "
                             << strsignal(WTERMSIG(res)));
  throw std::runtime_error(err::msg() << msg
                           << "\ncommand '" << command 
                           << "' terminated for an unknown reason");
}

void flock_throw(FILE* f, int op) {
     while(true) {
	  int res = flock(fileno(f), op);
	  if(res == 0)
	       return;
	  if(res != EINTR)
            throw std::runtime_error(err::msg() << "Failed to lock file: " << err::no);
     }
}

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
    return;
  int flags = fcntl(fd, F_GETFD);
  if(flags == -1)
    return;
  flags = fcntl(fd, F_SETFD, flags|FD_CLOEXEC);
  if(flags == -1)
    return;
}

int analyzeGap(struct arguments threadArg, FILE* resultFile, FILE* errFile);

cmdline_parse args; // GLOBAL: command line switches
std::string exeDir; // Location of executable

int main (int argc, char **argv)
{
  args.parse (argc, argv);

  struct arguments    threadArgs;
  FILE               *contigEndSeqFile;
  FILE               *meanAndStdevFile;
  FILE               *resultFile;
  FILE               *errFile;
  std::vector<FILE*>  readSeqFilesByDir;
  charb               tempBuffer(100);

  char *tempPtr = strrchr (argv[0], '/');
  if(tempPtr == NULL) {
    exeDir = std::string(".");
  } else {
    exeDir = std::string(argv[0], tempPtr - argv[0]);
  }

  resultFile = fopen(args.output_arg, "w");
  if(!resultFile)
    err::die(err::msg() << "Can't open output file '" << args.output_arg << "': " << err::no);
  close_on_exec(resultFile);

  errFile = fopen(args.error_out_arg, "a");
  if(!errFile)
    err::die(err::msg() << "Can't open error file: " << err::no);
  close_on_exec(errFile);

  contigEndSeqFile = fopen(args.contig_end_sequence_file_arg, "r");
  close_on_exec(contigEndSeqFile);
  if(!contigEndSeqFile)
    err::die(err::msg() << "Can't open '" << args.contig_end_sequence_file_arg << "': "
        << jellyfish::err::no);

  charb fn(200), line(100), tempLine(1000);

  for (int fileNum=0; 1; fileNum++) {
    sprintf (fn, "%s/readFile.%03d", args.dir_for_read_sequences_arg, fileNum);
    FILE *filetemp = fopen(fn, "r");
    if(!filetemp)
      break;
    close_on_exec(filetemp);
    readSeqFilesByDir.push_back (filetemp);
    fgets(line, 100, filetemp); // Gets past the first line ('>0')
  }


  // Here we make the output directory
  int res = mkdir(args.output_dir_arg, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if(res == -1 && errno != EEXIST)
    err::die(err::msg() << "Can't create output directory '" << args.output_dir_arg << err::no);


  meanAndStdevFile = NULL;
  if(args.mean_and_stdev_file_given) {
    meanAndStdevFile = fopen(args.mean_and_stdev_file_arg, "r");
    close_on_exec(meanAndStdevFile);
    if(!meanAndStdevFile)
      err::die(err::msg() << "Can't open file '" << args.mean_and_stdev_file_arg << "'"
          << err::no);
  }

  int		status	  = 0;
  long int	dirNum;
  bool		atEnd	  = false;
  unsigned int	fileIndex = 0;
  for (dirNum=0; ! atEnd; dirNum++) {
    if(dirNum >= args.num_threads_arg)
      wait(&status);

    threadArgs.clear();
    for (int i=0; i<args.num_joins_per_directory_arg; ++i) {
	 if (! fgets (tempLine, 1000, contigEndSeqFile))
	      break;
	 strcat (threadArgs.fauxReadFileDataStr, tempLine);
	 for (int j = 0; j < 3; ++j) {
	      fgets (tempLine, 1000, contigEndSeqFile);
	      strcat (threadArgs.fauxReadFileDataStr, tempLine);
	 }
    }
    //	       fgets_append (threadArgs.fauxReadFileDataStr, contigEndSeqFile); }

#if 0
    if (dirNum == 28)
	 asm ("int3;");
#endif
    atEnd = true;
    fileIndex=0;
    while (fileIndex < readSeqFilesByDir.size()) {
    int lineNum=0;
    while (fgets (line, 100, readSeqFilesByDir[fileIndex])) {
	 if (line[0] == '>') 
	      break;
	 if (lineNum % 2 == 0)
	      strcat (threadArgs.readFileDataStr, ">");
	 strcat (threadArgs.readFileDataStr, line);
	 ++lineNum;
    }
    ++fileIndex;
    }

    if(line[0] == '>') atEnd = false;
	 
    if (args.mean_and_stdev_file_given) {
      int dirNumTemp;
      int scanned = fscanf (meanAndStdevFile, "%d %d %d\n",
                            &dirNumTemp, &threadArgs.mean, &threadArgs.stdev);
      if(scanned != 3){
        std::cerr << "WARNING Bad mean and stdev for gap " << dirNumTemp << " " << threadArgs.mean << " " << threadArgs.stdev << " "  << jellyfish::err::no;
	threadArgs.mean  = 600;
	threadArgs.stdev = 200;
	}
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
      err::die(err::msg() << "Fork failed" << err::no);
    case 0:
	 analyzeGap(threadArgs, resultFile, errFile);
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
void do_analyzeGap(struct arguments& threadArg, const char* outDirName,
                   FILE* resultFile);
int analyzeGap(struct arguments threadArg, FILE* resultFile, FILE* errFile) {
  charb outDirName(100);

  sprintf (outDirName, "%s/gap%09lddir", args.output_dir_arg, threadArg.dirNum);
  int res = mkdir((char *)outDirName, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
  if(res == -1) {
    if(errno != EEXIST) { // OK if directory already exists
      fprintf(stderr, "Failed to create subdirectory '%s'\n", (const char*)outDirName);
      return -1;
    }
  }

  int ret_val = 0;
  try {
    do_analyzeGap(threadArg, outDirName, resultFile);
  } catch(std::runtime_error e) {
       flock_throw(errFile, LOCK_EX);
       fprintf(errFile, "Analyze gap failed for dir '%s': %s\n",
	       (const char*)outDirName, e.what());
       fflush(errFile);
       flock_throw(errFile, LOCK_UN);
       ret_val = -1;
  }

  if(!args.keep_directories_flag)
       rm_rf(outDirName);

  return ret_val;
}

void do_analyzeGap(struct arguments& threadArg, const char* outDirName,
                   FILE* resultFile) {
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
//  cmd << exeDir.c_str() << "/closeGaps.oneDirectory.perl"
  cmd << exeDir.c_str();
  if (args.jumping_read_joining_run_flag)
       cmd << "/closeGaps.oneDirectory.fromMinKmerLen.perl";
  else
       cmd << "/closeGaps.oneDirectory.perl";
  cmd << " 1>" << outDirName << "/out.err 2>&1"
      << " --dir-to-change-to " << outDirName
      << " --Celera-terminator-directory " << args.Celera_terminator_directory_arg
      << " --reads-file reads.fasta --output-directory outputDir"
      << " --max-kmer-len " << args.max_kmer_len_arg
      << " --min-kmer-len " << args.min_kmer_len_arg
      << " --maxnodes " << args.max_nodes_arg
      << " --mean-for-faux-inserts " << threadArg.mean
      << " --stdev-for-faux-inserts " << threadArg.stdev
      << " --num-stdevs-allowed " << args.num_stdevs_allowed_arg
      << " --join-aggressive " << args.join_aggressive_arg
      << " --use-all-kunitigs --noclean";

  system_throw(cmd.str().c_str());

  /* Now, if "passingKMer.txt" exists in outDirName, copy the files
     superReadSequences.fasta and
     readPlacementsInSuperReads.final.read.superRead.offset.ori.txt (after appropriate
     modifications), copy to the desired output directory from (e.g.)
     'outDirName'/work_localReadsFile_41_2 */
  if (args.num_joins_per_directory_arg > 1) {
       charb resultString, tempResultString;
       sprintf (tempFileName, "%s/passingReadsFile.txt", outDirName);
       infile = fopen_throw(tempFileName, "r", "Reading list of passing reads");
       if (!infile)
         throw std::runtime_error("Can't read the file containing the passing reads");
       while (fgets (tempResultString, 100, infile))
	    strcat (resultString, tempResultString);
       fclose (infile);
       flock_throw(resultFile, LOCK_EX);
       fputs((char *) resultString, resultFile);
       fflush (resultFile);
       flock_throw(resultFile, LOCK_UN);

       return;
  }


  sprintf (tempFileName, "%s/passingKMer.txt", outDirName);
  int passingKMerValue = 0;

  infile = fopen_throw(tempFileName, "r", "Reading passing k-mer size");
  if(!infile)
    throw std::runtime_error("Can't read passing k-mer value file");
  int scanned = fscanf(infile, "%d", &passingKMerValue);
  fclose(infile);
  if(scanned != 1)
    throw std::runtime_error(err::msg() << "Failed to read passing k-mer size"
                             << err::no);

  if(passingKMerValue < args.min_kmer_len_arg)
    throw std::runtime_error("Gap closing failed");

  charb superReadFastaString(1000);
  sprintf (tempFileName, "%s/work_localReadsFile_%d_2/superReadSequences.fasta",
           outDirName, passingKMerValue);
 infile = fopen_throw (tempFileName, "r", "Reading SuperRead sequences");
  fgets (line, 100, infile);
  sprintf (superReadFastaString, "%ld ", threadArg.dirNum);
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
    sprintf(readPlacementLine, " %s %ld %s %s", flds[0], threadArg.dirNum, flds[2], flds[3]);
    strcat(readPlacementLines,readPlacementLine);
  }
  fclose (infile);

  flock_throw(resultFile, LOCK_EX);
  fprintf(resultFile,"%s%s\n",(char*)superReadFastaString, (char*)readPlacementLines);
  fflush (resultFile);
  flock_throw(resultFile, LOCK_UN);

// #if 0
//   sprintf (cmd, "cp %s/work_localReadsFile_%d_2/superReadSequences.fasta %s/superReadSequences.%09d.fasta", outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
//   system (cmd);
//   sprintf (cmd, "cp %s/work_localReadsFile_%d_2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt %s/readPlacementsInSuperReads.final.read.superRead.offset.ori.%09d.txt", outDirName, passingKMerValue, args.output_dir_arg, threadArg.dirNum);
//   system (cmd);
// #endif

  return;
}
