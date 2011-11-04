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

#include <jellyfish/err.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/invertible_hash_array.hpp>
#include <jellyfish/allocators_mmap.hpp>
#include <stdint.h>

#include <iostream>

#include "diskBasedUnitigger.h"
#define KUNITIG_FILE "/genome8/raid/tri/kUnitigStudy/arg_ant/afterAlekseyAndMikeRedundentKill/guillaumeKUnitigsAtLeast32bases_all.fasta"
#define READ_DATA_FILE "/genome8/raid/tri/testDirForReadPlacementRoutines/brucellaData/brucella.pass5reads.fasta"

unsigned int  *kunitigNumber, *kunitigOffset;
unsigned long *kUnitigLengths;
uint64_t       hash_size;
int            lastKUnitigNumber;
unsigned char *kmerOriInKunitig;
char           line[1000000], kUnitigBases[1000000], kUnitigBasesRevComp[1000000];
uint64_t       kUnitigKmers[1000000], kUnitigRevCompKmers[1000000], mask;
int            readDataFile;
char           hdrLine[10000];
int            readNumber;
int            numFilenames, longOutput;
char          *filenames[1000];
char           charToBaseValue[256];
uint64_t       charToRevCompBaseValue[256];
int            mer_len;
int            nb_threads = 5;
const char    *prefix = 0;

/* Our types.
 */
//typedef jellyfish::invertible_hash::array<uint64_t,atomic::gcc<uint64_t>,allocators::mmap> inv_hash_storage_t;
typedef inv_hash_storage_t::iterator iterator_t;

void initializeValues (void);
void getMatchesForRead (const char *readBasesBegin, const char *readBasesEnd, 
                        const char *readName, inv_hash_storage_t *hashPtr, FILE *out);


class ProcessReads : public thread_exec {
  jellyfish::parse_read  read_parser;
  inv_hash_storage_t     *hash;
  const char             *prefix;

public:
  ProcessReads(int argc, char *argv[], inv_hash_storage_t *h, const char *pref) : 
    read_parser(argc, argv, 100), hash(h), prefix(pref) {}

  virtual void start(int id) {
    jellyfish::parse_read::thread  read_stream(read_parser.new_thread());
    jellyfish::parse_read::read_t *read;
    char                           readBases[100000];
    char                           header[3000];
    char                           out_file[3000];
    FILE                          *out;
    if(!prefix) {
      out = stdout;
    } else {
      sprintf(out_file, "%s_%d", prefix, id);
      out = fopen(out_file, "w");
      if(!out)
        die << "Can't open output file '" << out_file << "'" << err::no;
    }

    while((read = read_stream.next_read())) {
      strncpy(header, read->header, read->hlen);
      header[read->hlen] = '\0';
      strtok(header, " ");
      char *optr = readBases;
      for(const char *iptr = read->seq_s; iptr < read->seq_e; iptr++) {
        if(!isspace(*iptr)) {
          *optr = *iptr;
          optr++;
        }
      }
      getMatchesForRead(readBases, optr, header, hash, out);
    }
  }
};

int main(int argc, char *argv[])
{
  FILE *infile;
  int i, kUnitigNumber;
  char *kUnitigFilename;
  uint64_t searchValue;
  uint64_t val, id;
  char *cptr;
  //  char readName[1000];
  char *numKUnitigsFile;
  int kUnitigLength;
  unsigned char ori;

  /* Database file name is first argument. The k-unitig filename
     is the second; The file with the number of k-unitigs is third;
     after that follow the read file(s)
  */
  longOutput = 0;
  filenames[1] = (char *) KUNITIG_FILE;
  filenames[3] = (char *) READ_DATA_FILE;
  numFilenames = 0;
  for (i=1; i<argc; i++) {
    if (strcmp(argv[i], "-l") == 0) {
      longOutput = 1;
      continue; }
    if(strcmp(argv[i], "-p") == 0) {
      ++i;
      prefix = argv[i];
      continue;
    }
    if(strcmp(argv[i], "-t") == 0) {
      ++i;
      nb_threads = atoi (argv[i]);
      continue;
    }
    filenames[numFilenames] = argv[i];
    ++numFilenames; }

  numKUnitigsFile = filenames[2];
  kUnitigFilename = filenames[1];
  // Get input data reads directly from filenames
  // if (numFilenames > 3)
  //   readDataFile = filenames[3];

  /* Map the hash in memory. Yeah, it should/could be shorter.
   */
  mapped_file dbf(filenames[0]);
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
     
  initializeValues ();

  // create our mask
  mask = 0;
  for (i=0; i<mer_len; i++)
    mask = (mask << 2) + 0x3;
  // Find out the last kUnitig number
  infile = fopen (numKUnitigsFile, "r");
  if(!infile)
    die << "Failed to open file '" << numKUnitigsFile << "'" << err::no;
  int fields_read = fscanf (infile, "%d\n", &lastKUnitigNumber);
  if(fields_read != 1)
    die << "Failed to read the last k-unitig number from file '"
        << numKUnitigsFile << "'" << err::no;
  fclose (infile);
  fprintf (stderr, "The largest kUnitigNumber was %d\n", lastKUnitigNumber);
  mallocOrDie (kUnitigLengths, (lastKUnitigNumber+2), uint64_t);
     
  fprintf (stderr, "Opening file %s...\n", kUnitigFilename);
  infile = fopen (kUnitigFilename, "r");
  if(!fgets (line, sizeof(line), infile)) // This is a header line
    die << "Failed to read header line from file '"
        << kUnitigFilename << "'" << err::no;
  fields_read = sscanf (line, ">%d length:%d", &kUnitigNumber, &kUnitigLength);
  if(fields_read != 2)
    die << "Header of file '" << kUnitigFilename
        << "' does not match pattern '>UnitigiNumber length:UnitigLength'";
  kUnitigLengths[kUnitigNumber] = kUnitigLength;
  if (kUnitigNumber % 100 == 0)
    fprintf (stderr, "\rkUnitigNumber = %d", kUnitigNumber);
  //printf ("kUnitigNumber = %d; kUnitigLength = %d\n", kUnitigNumber, kUnitigLength);
  cptr = kUnitigBases;
  while (1) {
    int atEof = 0;
    if (! fgets(line, sizeof(line), infile)) {
      atEof = 1;
      goto processTheKUnitig; }
    if (line[0] == '>')
      goto processTheKUnitig;
    strcpy (cptr, line);
    while (! isspace(*cptr)) {
      switch (*cptr) {
      case 'a': case 'A': *cptr = 0; break;
      case 'c': case 'C': *cptr = 1; break;
      case 'g': case 'G': *cptr = 2; break;
      case 't': case 'T': *cptr = 3; break; }
      ++cptr; }
    continue;
  processTheKUnitig:
    for (int j=0, k=kUnitigLength-1; j<kUnitigLength; j++, k--)
      kUnitigBasesRevComp[k] = 3 - kUnitigBases[j];
    // The kUnitigRevCompKmers are where they are placed in the
    // reverse complement; we will need to adjust at the end to
    // place them in relation to the forward direction of the k-unitig
    kUnitigKmers[0] = kUnitigRevCompKmers[0] = 0;
    for (int j=0; j<mer_len; j++) {
      kUnitigKmers[0] = (kUnitigKmers[0] << 2) + kUnitigBases[j];
      //	    printf ("j=%d, kmer=%x; ", j, kUnitigKmers[0]);
      kUnitigRevCompKmers[0] = (kUnitigRevCompKmers[0] << 2) + kUnitigBasesRevComp[j]; }
    //     printf ("\n");
    for (int j=mer_len, k=1; j<kUnitigLength; j++, k++) {
      kUnitigKmers[k] = (kUnitigKmers[k-1] << 2) + kUnitigBases[j];
      kUnitigKmers[k] &= mask;
      kUnitigRevCompKmers[k] = (kUnitigRevCompKmers[k-1] << 2) + kUnitigBasesRevComp[j];
      kUnitigRevCompKmers[k] &= mask; }
    if (kUnitigNumber % 100 == 0)
      fprintf (stderr, "\rkUnitigNumber = %d", kUnitigNumber);
    for (int j=0, k=kUnitigLength-mer_len; k>=0; j++, k--) {
      if (kUnitigKmers[j] < kUnitigRevCompKmers[k]) {
        searchValue = kUnitigKmers[j];
        ori = 0; }
      else {
        searchValue = kUnitigRevCompKmers[k];
        ori = 1; }
      //	    fprintf (stdout, "offset = %d, forward mer value = %llx, backwards mer value = %llx\n", j, kUnitigKmers[j], kUnitigRevCompKmers[k]);
      // Here you have to search for the index, then install the
      // information into your arrays at the index
      if (! hash->get_val(searchValue, id, val, false)) {
        fprintf (stderr, "mer %lx was not found. Bye!\n", searchValue);
        exit (1); // WAS OUT FOR DEBUGGING
      }
      //	    fprintf (stdout, "id = %lld, val = %lld\n", id, val);
      kunitigNumber[id] = kUnitigNumber;
      kunitigOffset[id] = j;
      kmerOriInKunitig[id] = ori;
    }
    //     break; // For debugging
    if (atEof)
      break;
    sscanf (line, ">%d length:%d", &kUnitigNumber, &kUnitigLength);
    kUnitigLengths[kUnitigNumber] = kUnitigLength;
    //     printf ("kUnitigNumber = %d; kUnitigLength = %d\n", kUnitigNumber, kUnitigLength);
    cptr = kUnitigBases;
  }
  //exit (0); // For debugging

  ProcessReads process_reads(numFilenames - 3, filenames + 3, hash, prefix);
  process_reads.exec_join(nb_threads);
     
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
                        FILE *out)
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
    fprintf (out, "readNumber = %d, readLength = %d\n", readNumber, readLength);
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
  }
}

