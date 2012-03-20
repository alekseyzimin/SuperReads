/* Arguments are as follows:
   1) The name of the database output by jellyfish (using raw data)
   2) The name of the file with the k-unitigs
   3) The name of the file which reports how many k-unitigs there are
   4) The name of the file containing the reads

   4/12/11: The output has now been reduced; to get the longer (former output),
   use the '-l' flag.
   5/11/11: Made parallel
   5/11/11: Use the '-p' flag to give an output prefix.
   5/13/11: Use the '-t' flag to specify how many threads should run.
*/

#ifndef typeof
#define typeof __typeof__
#endif

#include <src/mer_dna.hpp>
#include <src/read_parser.hpp>

// #include <jellyfish/err.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/invertible_hash_array.hpp>
#include <jellyfish/allocators_mmap.hpp>
#include <stdint.h>
#include <cstdio>

#include <iostream>
#include <fstream>

#include <src/diskBasedUnitigger.h>
#define KUNITIG_FILE "/genome8/raid/tri/kUnitigStudy/arg_ant/afterAlekseyAndMikeRedundentKill/guillaumeKUnitigsAtLeast32bases_all.fasta"
#define READ_DATA_FILE "/genome8/raid/tri/testDirForReadPlacementRoutines/brucellaData/brucella.pass5reads.fasta"

#include <charb.hpp>
#include <gzip_stream.hpp>
#include <jflib/multiplexed_io.hpp>

#include <src2/findMatchesBetweenKUnitigsAndReads.hpp>

unsigned int  *kunitigNumber, *kunitigOffset;
unsigned long *kUnitigLengths;
uint64_t       hash_size;
int            lastKUnitigNumber;
unsigned char *kmerOriInKunitig;
charb          line(1000);
uint64_t       mask;
int            readNumber; // TODO: Remove: it is not really used
int            longOutput;
char           charToBaseValue[256];
uint64_t       charToRevCompBaseValue[256];
int            mer_len;

/* Our types.
 */
//typedef jellyfish::invertible_hash::array<uint64_t,atomic::gcc<uint64_t>,allocators::mmap> inv_hash_storage_t;
typedef inv_hash_storage_t::iterator iterator_t;
class header_output;

void initializeValues (void);
void getMatchesForRead (const char *readBasesBegin, const char *readBasesEnd, 
                        const char *readName, inv_hash_storage_t *hashPtr, header_output& out);

class ReadKunitigs : public thread_exec {
  read_parser         unitig_parser;
  inv_hash_storage_t* hash;
  
public:
  ReadKunitigs(const char* kUnitigFile, inv_hash_storage_t* hashPtr, int nb_threads) :
    unitig_parser(kUnitigFile, nb_threads), hash(hashPtr) { }

  virtual void start(int th_id) {
    read_parser::stream unitig_stream(unitig_parser);
    mer_dna fwd_mer(mer_len), rev_mer(mer_len);
    for( ; unitig_stream; ++unitig_stream) {
      // Parse header
      int kUnitigNumber, kUnitigLength;
      int fields_read = sscanf(unitig_stream->header, ">%d length:%d", &kUnitigNumber, &kUnitigLength);
      kUnitigLengths[kUnitigNumber] = kUnitigLength;
      if(fields_read != 2)
        die << "Fasta header of file does not match pattern '>UnitigiNumber length:UnitigLength'\n"
            << unitig_stream->header << "\n";
      if(unitig_stream->sequence.size() < (size_t)mer_len)
        continue;
      char* cptr = unitig_stream->sequence;
      // Prime the k-mers
      for(int i = 0; i < mer_len - 1; ++i, ++cptr) {
        fwd_mer.shift_left(*cptr);
        rev_mer.shift_right(0x3 - (fwd_mer[0] & 0x3));
      }

      for( ; cptr < unitig_stream->sequence.end(); ++cptr) {
        fwd_mer.shift_left(*cptr);
        rev_mer.shift_right(0x3 - (fwd_mer[0] & 0x3));
         
        uint64_t searchValue;
        unsigned char ori;
        if(fwd_mer[0] < rev_mer[0]) {
          searchValue = fwd_mer[0];
          ori = 0;
        } else {
          searchValue = rev_mer[0];
          ori = 1;
        }
        uint64_t id = 0, val;
        if(!hash->get_val(searchValue, id, val, false)) {
          fprintf (stderr, "mer %lx was not found. Bye!\n", searchValue);
          //		    exit (1); // WAS OUT FOR DEBUGGING
        }
        kunitigNumber[id]             = kUnitigNumber;
        kunitigOffset[id]             = cptr - unitig_stream->sequence.begin() - (mer_len - 1);
        kmerOriInKunitig[id]          = ori;
      }
    }
  }
};

/** Keep track whether or not a header has been written. On the first
    actual output (operator<<), the header is prepended. set_header()
    must be called to set/create the header. After we are done with a
    read, flush() is called, which adds a new line if necessary and
    endr (end of record) if asked to.
 */
class header_output {
  jflib::omstream& out_;
  charb            header_;
  bool             is_written_; // whether header has been written already
public:
  header_output(jflib::omstream& out) : out_(out), is_written_(false) { }

  void set_header(const char* name, size_t length) {
    sprintf(header_, "%s %ld", name, length);
    is_written_ = false;
  }
  void flush(bool end_of_record) {
    if(is_written_)
      out_ << "\n";
    if(end_of_record)
      out_ << jflib::endr;
  }

  template<typename T>
  jflib::omstream& operator<<(T& x) {
    if(!is_written_) {
      out_ << header_;
      is_written_ = true;
    }
    out_ << x;
    return out_;
  }
};

class ProcessReads : public thread_exec {
  read_parser           parser;
  inv_hash_storage_t   *hash;
  jflib::o_multiplexer  multiplexer;

public:
  template<typename Iterator>
  ProcessReads(Iterator file_start, Iterator file_end, 
               inv_hash_storage_t *h, std::ostream& out_stream, int nb_threads) : 
    parser(file_start, file_end, nb_threads), 
    hash(h),
    multiplexer(&out_stream, 3 * nb_threads, 4096)
  {}

     virtual void start(int id) {
       read_parser::stream read_stream(parser);
       jflib::omstream     out(multiplexer);
       header_output       hout(out);
       char                read_prefix[3], prev_read_prefix[3];
       uint64_t            read_id = 0, prev_read_id = 0;
       memset(read_prefix, '\0', sizeof(read_prefix));
       
       for( ; read_stream; ++read_stream) {
         memcpy(prev_read_prefix, read_prefix, sizeof(read_prefix));
         prev_read_id = read_id;
         
         strtok(read_stream->header, " ");
         sscanf(read_stream->header, ">%2s%ld", read_prefix, &read_id);
         
         // Start new line & print read header if new read
         bool is_new_read = strcmp(read_prefix, prev_read_prefix) || read_id != prev_read_id;
         if(is_new_read) {
           // Keep mated read together in output.
           bool end_of_record =
             (read_id % 2 == 0) || (read_id != (prev_read_id + 1)) || strcmp(read_prefix, prev_read_prefix);
           hout.flush(end_of_record);
           hout.set_header(read_stream->header + 1, read_stream->sequence.size());
         }

         getMatchesForRead(read_stream->sequence, read_stream->sequence.end(),
                           read_stream->header + 1, hash, hout);
       }
       hout.flush(true);
     }
};

int main(int argc, char *argv[])
{
     FILE                               *infile;
     findMatchesBetweenKUnitigsAndReads  args(argc, argv);

     longOutput = args.long_flag;

     mapped_file dbf(args.jellyfishdb_arg);
     if(memcmp(dbf.base(), "JFRHSHDN", 8))
	  die << "Invalid database format, expected 'JFRHSHDN'";
     dbf.random().load();
     raw_inv_hash_query_t qhash(dbf);
     //  inv_hash_storage_t hash(dbf.base() + 8, dbf.length() - 8);
     inv_hash_storage_t *hash = qhash.get_ary();
     
     /* Get some info about the hash.
      */
     //  unsigned int mer_len = hash->get_key_len() / 2;
     mer_len = hash->get_key_len() / 2;
     std::cerr << "mer length: " << mer_len << "\n";
     std::cerr << "hash size: " << hash->get_size() << "\n";
     
     hash_size = hash->get_size();
     // Allocate the space for the other data
     mallocOrDie (kunitigNumber, hash_size, unsigned int);
     mallocOrDie (kunitigOffset, hash_size, unsigned int);
     mallocOrDie (kmerOriInKunitig, hash_size, unsigned char);

     // Open output file
     gzipstream out(args.output_arg);
     
     initializeValues ();

     // create our mask
     mask = ((uint64_t)1 << (mer_len * 2)) - 1;

     // Find out the last kUnitig number
     infile = fopen (args.numKUnitigsFile_arg, "r");
     if(!infile)
	  die << "Failed to open file '" << args.numKUnitigsFile_arg << "'" << err::no;
     int fields_read = fscanf (infile, "%d\n", &lastKUnitigNumber);
     if(fields_read != 1)
	  die << "Failed to read the last k-unitig number from file '"
	      << args.numKUnitigsFile_arg << "'" << err::no;
     fclose (infile);
     fprintf (stderr, "The largest kUnitigNumber was %d\n", lastKUnitigNumber);
     mallocOrDie (kUnitigLengths, (lastKUnitigNumber+2), uint64_t);

     { // Read k-unitigs k-mers and positions
       ReadKunitigs KUnitigReader(args.kUnitigFile_arg, hash, args.threads_arg);
       KUnitigReader.exec_join(args.threads_arg);
     }
     

     ProcessReads process_reads(args.readFiles_arg.begin(), args.readFiles_arg.end(),
                                hash, out, args.threads_arg);
     process_reads.exec_join(args.threads_arg);
     
     // Now working with the reads
     // fprintf (stderr, "Opening the read file: %s...\n", readDataFile);
     // readNumber = 0;
     // infile = fopen (readDataFile, "r");
     // fgets (hdrLine, 10000, infile); // This is a header line
     // //  sscanf (line, ">%d length:%d", &kUnitigNumber, &kUnitigLength);
     // sscanf (hdrLine, ">%s", readName);
     // //  fprintf (stdout, "readNumber = %d", readNumber);
     // //       printf ("kUnitigNumber = %d; kUnitigLength = %d\n", kUnitigNumber, kUnitigLength);
     // cptr = readBases;
     // while (1) {
     //   int atEof = 0;
     //   if (! fgets(cptr, 10000, infile)) {
     //     atEof = 1;
     //     goto processTheRead; }
     //   if (cptr[0] == '>')
     //     goto processTheRead;
     //   while (! isspace(*cptr)) {
     //     switch (*cptr) {
     //     case 'a': case 'A':
     //     case 'c': case 'C':
     //     case 'g': case 'G':
     //     case 't': case 'T':
     //       break;
     //     default:
     //       *cptr = 'N';
     //       break;
     //     }
     //     ++cptr; }
     //   continue;
     // processTheRead:
     //   getMatchesForRead (readBases, cptr, readName, hash);
     //   if (atEof)
     //     break;
     //   strcpy (hdrLine, cptr);
     //   sscanf (hdrLine, ">%s", readName);
     //   ++readNumber;
     //   cptr = readBases;
     // }
     
     return (0);
}

void initializeValues (void)
{
     memset (charToBaseValue, 0, 256);
     charToBaseValue[(int)'a'] = 0;
     charToBaseValue[(int)'A'] = 0;
     charToBaseValue[(int)'c'] = 1;
     charToBaseValue[(int)'C'] = 1;
     charToBaseValue[(int)'g'] = 2;
     charToBaseValue[(int)'G'] = 2;
     charToBaseValue[(int)'t'] = 3;
     charToBaseValue[(int)'T'] = 3;
     memset (charToRevCompBaseValue, 0, 256 * sizeof (uint64_t));
     charToRevCompBaseValue[(int)'a'] = 3ULL << (2 * (mer_len-1));
     charToRevCompBaseValue[(int)'A'] = 3ULL << (2 * (mer_len-1));
     charToRevCompBaseValue[(int)'c'] = 2ULL << (2 * (mer_len-1));
     charToRevCompBaseValue[(int)'C'] = 2ULL << (2 * (mer_len-1));
     charToRevCompBaseValue[(int)'g'] = 1ULL << (2 * (mer_len-1));
     charToRevCompBaseValue[(int)'G'] = 1ULL << (2 * (mer_len-1));
     charToRevCompBaseValue[(int)'t'] = 0ULL << (2 * (mer_len-1));
     charToRevCompBaseValue[(int)'T'] = 0ULL << (2 * (mer_len-1));
}

void getMatchesForRead (const char *readBasesBegin, const char *readBasesEnd, 
                        const char *readName, inv_hash_storage_t *hashPtr,
                        header_output& out)
{
     int readLength;
     uint64_t readKmer, readKmerTmp, readRevCompKmer, readRevCompKmerTmp;
     uint64_t searchValue;
     unsigned char kUnitigOri, netOri, netOriHold, ori;
     int nextWithN;
     int jtmp, ktmp, tempIndex;
     int kUnitigNumber, kUnitigOffset;
     int isFirstRecord;
     int kUnitigOffsetOfFirstBaseInRead, kUnitigOffsetOfLastBaseInRead;
     int kUnitigOffsetOfFirstBaseInReadHold=0, kUnitigNumberHold=0;
     int intervalBegin=0, intervalEnd=0;
     int ahg=0, bhg=0;
     int firstKUnitigOffsetForNucmer=0, lastKUnitigOffsetForNucmer=0;
     int firstReadOffsetForNucmer=0, lastReadOffsetForNucmer=0;
     uint64_t val, id;

     readLength = readBasesEnd - readBasesBegin;
     if (longOutput)
         out << "readNumber = " << readNumber 
             << " readLength = " << readLength << "\n";
       //fprintf (out, "readNumber = %d, readLength = %d\n", readNumber, readLength);
     readKmer = readRevCompKmer = 0;
     if (readLength >= mer_len)
	  for (int j=0; j<mer_len; j++) {
	       readKmer = (readKmer << 2) + charToBaseValue[(int)readBasesBegin[j]];
	       readRevCompKmer = (readRevCompKmer >> 2) + charToRevCompBaseValue[(int)readBasesBegin[j]]; }
     netOriHold = 2;
     nextWithN = readLength;
     for (int j=0; j<readLength; j++)
	  if (readBasesBegin[j] == 'N') {
	       nextWithN = j;
	       break; }
     for (int j=0, k=readLength-mer_len; k>=0; j++, k--) {
	  // Re-setting nextWithN if necessary
	  if (nextWithN < j)
	       for (jtmp=j; jtmp<readLength; jtmp++)
		    if (readBasesBegin[jtmp] == 'N') {
			 nextWithN = jtmp;
			 break; }
	  if (nextWithN < j)
	       nextWithN = readLength;
	  if (j>0) {
	       readKmer = (readKmer << 2) + charToBaseValue[(int)readBasesBegin[j+mer_len-1]];
	       readRevCompKmer = (readRevCompKmer >> 2) + charToRevCompBaseValue[(int)readBasesBegin[j+mer_len-1]];
	       readKmer &= mask;
	       readRevCompKmer &= mask; }
	  if (nextWithN-j < mer_len)
	       continue;
	  if (readKmer < readRevCompKmer) {
	       searchValue = readKmer;
	       ori = 0; }
	  else {
	       searchValue = readRevCompKmer;
	       ori = 1; }
	  if (! hashPtr->get_val(searchValue, id, val, false)) {
	       continue;
	  }
	  // It appears sometimes it comes back that the k-mer is
	  // in the data but the kUnitigNumber and offset are 0
	  // We must eliminate these records.
	  // (Checking below)
	  kUnitigNumber = kunitigNumber[id];
	  kUnitigOffset = kunitigOffset[id];
	  kUnitigOri = kmerOriInKunitig[id];
	  if ((kUnitigNumber == 0) && (kUnitigOffset == 0))
	       continue;
	  // now work with ori, j, and readNumber
	  isFirstRecord = 0;
	  if (kUnitigOri == ori)
	       netOri = 0;
	  else
	       netOri = 1;
	  // in the following kUnitigOffsetOfFirstBaseInRead
	  if (netOri == 0) // The read is 'F' in kUnitig
	       kUnitigOffsetOfFirstBaseInRead = kUnitigOffset - j;
	  else
	       kUnitigOffsetOfFirstBaseInRead = kUnitigOffset + j + (mer_len);
	  if ((kUnitigOffsetOfFirstBaseInReadHold != kUnitigOffsetOfFirstBaseInRead) ||
	      (kUnitigNumber != kUnitigNumberHold) ||
	      (netOri != netOriHold) ||
	      (j > intervalEnd)) {
	       if (netOriHold < 2) {
#if 0
		    if (longOutput)
			 fprintf (out, "Portion of read covered by the last k-unitig match ");
		    fprintf (out, "= (%d, %d)", intervalBegin, intervalEnd);
		    if (longOutput)
			 fprintf (out, "\nmyNucmerLine (ahg,bhg) =");
		    fprintf (out, " (%d,%d): %d %d", ahg, bhg, firstKUnitigOffsetForNucmer, lastKUnitigOffsetForNucmer);
		    if (netOriHold == 0)
			 fprintf (out, " %d %d", firstReadOffsetForNucmer, lastReadOffsetForNucmer);
		    else
			 fprintf (out, " %d %d", lastReadOffsetForNucmer, firstReadOffsetForNucmer);
		    fprintf (out, " %lu %d %d %s", kUnitigLengths[kUnitigNumberHold], readLength, kUnitigNumberHold, readName);
		    fprintf (out, "\n");
#else
                    out << " " << kUnitigNumberHold
                        << " " << ahg
                        << " " << ((netOriHold == 0) ? 'F' : 'R');
                    //		    fprintf (out, "%c %d %d %d %s\n", (netOriHold == 0) ? 'F' : 'R', ahg, readLength, kUnitigNumberHold, readName);
#endif
	       }
	       kUnitigOffsetOfFirstBaseInReadHold = kUnitigOffsetOfFirstBaseInRead;
	       kUnitigNumberHold = kUnitigNumber;
	       netOriHold = netOri;
	       intervalBegin = j; 
	       isFirstRecord = 1; }
	  intervalEnd = j + mer_len;
	  if (netOri == 0)
	       kUnitigOffsetOfLastBaseInRead = kUnitigOffsetOfFirstBaseInRead + readLength;
	  else
	       kUnitigOffsetOfLastBaseInRead = kUnitigOffsetOfFirstBaseInRead - readLength;
	  // Do Celera-style end calculations
	  // Use the k-unitig as read1 and the read as read2
	  if (netOri == 0) {
	       ahg = kUnitigOffsetOfFirstBaseInRead;
	       bhg = kUnitigOffsetOfLastBaseInRead - kUnitigLengths[kUnitigNumber]; }
	  else {
	       ahg = kUnitigOffsetOfLastBaseInRead;
	       bhg = kUnitigOffsetOfFirstBaseInRead - kUnitigLengths[kUnitigNumber]; }
	  ////////// BEGIN NUCMER END CALCULATIONS
	  if (ahg >= 0) {
	       firstKUnitigOffsetForNucmer = ahg + 1;
	       if (netOri == 0)
		    firstReadOffsetForNucmer = 1;
	       else
		    lastReadOffsetForNucmer = readLength;
	  }
	  else {
	       firstKUnitigOffsetForNucmer = 1;
	       if (netOri == 0)
		    firstReadOffsetForNucmer = 1 - ahg;
	       else
		    lastReadOffsetForNucmer = readLength + ahg;
	  }
	  if (bhg  >= 0) {
	       lastKUnitigOffsetForNucmer = kUnitigLengths[kUnitigNumber];
	       if (netOri == 0)
		    lastReadOffsetForNucmer = readLength - bhg;
	       else
		    firstReadOffsetForNucmer = 1 + bhg;
	  }
	  else {
	       lastKUnitigOffsetForNucmer = kUnitigLengths[kUnitigNumber] + bhg;
	       if (netOri == 0)
		    lastReadOffsetForNucmer = readLength;
	       else
		    firstReadOffsetForNucmer = 1;
	  }
	  ////////// END NUCMER END CALCULATIONS
	  // Do we match at the end of the possible interval? (can jump)
	  if (isFirstRecord) {
	       jtmp = lastReadOffsetForNucmer - mer_len;
	       ktmp = readLength - jtmp - mer_len;
	       // If the desired k-mer has an 'N', don't bother
	       if (nextWithN < jtmp) {
		    for (tempIndex=jtmp; tempIndex<jtmp+mer_len; tempIndex++)
			 if (readBasesBegin[tempIndex] == 'N')
			      break; }
	       else
		    tempIndex = nextWithN;
	       if (tempIndex-jtmp < mer_len) // It had an 'N'
		    continue;
	       readKmerTmp = readKmer;
	       readRevCompKmerTmp = readRevCompKmer;
	       tempIndex = j + mer_len - jtmp;
	       if (tempIndex < 0)
		    tempIndex = 0;
	       for (tempIndex=0; tempIndex<mer_len; tempIndex++) {
		    readKmerTmp = (readKmerTmp << 2) + charToBaseValue[(int)readBasesBegin[jtmp+tempIndex]];
		    readRevCompKmerTmp = (readRevCompKmerTmp >> 2) + charToRevCompBaseValue[(int)readBasesBegin[jtmp+tempIndex]]; }
	 
	       readKmerTmp &= mask;
	       readRevCompKmerTmp &= mask;
	 
	       if (readKmerTmp < readRevCompKmerTmp) {
		    searchValue = readKmerTmp;
		    ori = 0; }
	       else {
		    searchValue = readRevCompKmerTmp;
		    ori = 1; }
	       if (! hashPtr->get_val(searchValue, id, val, false)) {
		    continue;
	       }
	       // It appears sometimes it comes back that the k-mer is
	       // in the data but the kUnitigNumber and offset are 0
	       // We must eliminate these records.
	       // (Checking below)
	       kUnitigNumber = kunitigNumber[id];
	       kUnitigOffset = kunitigOffset[id];
	       kUnitigOri = kmerOriInKunitig[id];
	       if ((kUnitigNumber == 0) && (kUnitigOffset == 0))
		    continue;
	       // now work with ori, j, and readNumber
	       if (kUnitigOri == ori)
		    netOri = 0;
	       else
		    netOri = 1;
	       // in the following kUnitigOffsetOfFirstBaseInRead
	       if (netOri == 0) // The read is 'F' in kUnitig
		    kUnitigOffsetOfFirstBaseInRead = kUnitigOffset - jtmp;
	       else
		    kUnitigOffsetOfFirstBaseInRead = kUnitigOffset + jtmp + (mer_len);
	       if ((kUnitigOffsetOfFirstBaseInReadHold == kUnitigOffsetOfFirstBaseInRead) &&
		   (kUnitigNumber == kUnitigNumberHold) &&
		   (netOri == netOriHold)) {
		    j = jtmp;
		    k = ktmp; 
		    readKmer = readKmerTmp;
		    readRevCompKmer = readRevCompKmerTmp;
		    intervalEnd = j + mer_len; }
	  }
     }
     //	    break; // For debugging
     if (netOriHold < 2) {
#if 0
	  if (longOutput)
	       fprintf (out, "Portion of read covered by the last k-unitig match ");
	  fprintf (out, "= (%d, %d)", intervalBegin, intervalEnd);
	  if (longOutput)
	       fprintf (out, "\nmyNucmerLine (ahg,bhg) =");
	  fprintf (out, " (%d,%d): %d %d", ahg, bhg, firstKUnitigOffsetForNucmer, lastKUnitigOffsetForNucmer);
	  if (netOriHold == 0)
	       fprintf (out, " %d %d", firstReadOffsetForNucmer, lastReadOffsetForNucmer);
	  else
	       fprintf (out, " %d %d", lastReadOffsetForNucmer, firstReadOffsetForNucmer);
	  fprintf (out, " %lu %d %d %s", kUnitigLengths[kUnitigNumberHold], readLength, kUnitigNumberHold, readName);
	  fprintf (out, "\n");
#else
          //	  fprintf (out, "%c %d %d %d %s\n", (netOriHold == 0) ? 'F' : 'R', ahg, readLength, kUnitigNumberHold, readName);
          out << " " << kUnitigNumberHold
              << " " << ahg
              << " " << ((netOriHold == 0) ? 'F' : 'R');
#endif
     }
}

