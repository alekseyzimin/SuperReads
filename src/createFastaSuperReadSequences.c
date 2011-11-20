/*
I)   Allocate space for the k-unitigs
        A) Use the length of guillaumeKUnitigsAtLeast32bases_all.fasta
II)  Read in the file with the k-unitigs
        A) guillaumeKUnitigsAtLeast32bases_all.fasta
III) Allocate space for the pointers to the k-unitigs
        A) Find the largest number assigned to a k-unitig
              1) Use the files guillaumeKUnitigsAtLeast32bases_*.fa
IV)  Process the file (in memory) with the k-unitigs
        A) Set the pointers to the k-unitigs
        B) Eliminate spaces from the k-unitigs
V)   For each super-read
        A) Create the header
              1) Use the line with the order of the k-unitigs, their oris, and
                   their overlaps
                    a) No space between the k-unitig number and its orientation
                    b) Underscores between all other fields
        B) Generate the sequence
              1) Create the reverse-complement of a k-unitig if necessary
                    a) Must generate space if necessary ahead of time
              2) Can use strcat and strcpy
              3) Keep track of the output length
        C) Output the fastq format for the read
        D) Output quals 'a' for the length of the read
VI)  Make sure you allow for an appropriate set of params
        A) Working directory (where to find k-unitig sequence files)
        B) File listing the super-reads (with path)
              1) allData/superReadGroups.onePerLine.withReadInfoIncluded.txt
        C) -seqdiffmax # allows one to specify the number of differences one
	     allows between the sequences of overlapping k-unitigs in a
	     super-read; otherwise it fails (default is 0).
	D) -nosequence just outputs the names of the passing super-reads
	     instead of the name and sequence in fasta format
	E) -error-filename filename sends the error output to the specified
	     file
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>

#define NUM_KUNITIGS_FILENAME "numKUnitigs.txt"
#define KUNITIG_FILE_COMPLETE "guillaumeKUnitigsAtLeast32bases_all.fasta"
#define DEFAULT_SUPER_READ_LIST_FILE "superReadGroups.onePerLine.withReadInfoIncluded.txt"
#define AFTER_NEWLINE 1
#define AFTER_HEADER_NEWLINE 2
#define IN_SEQUENCE 3
#define MAX_READ_LEN 100000000

char *kUnitigSpace, **kUnitigSeq;
int *kUnitigLengths;
char *line, *reverseComplementSpace, *outputSeqSpace;

void generateReverseComplement (char *seq, int seqLen);
FILE *Fopen (const char *fn, const char *mode);

#define processTheChar if (! isspace (kUnitigSpace[i64])) { \
		           if (i64 != j64) \
			       kUnitigSpace[j64] = kUnitigSpace[i64]; \
		               ++j64; } \
	               state = IN_SEQUENCE; \
	               continue;

#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

int main (int argc, char **argv)
{
     char *workingDir, superReadListFile[500], fname[500], kUnitigFilename[500];
     struct stat statbuf;
     uint64_t kUnitigSeqFileSize, fsize;
     uint64_t i64, j64=0;
     char errorFilename[500];
     FILE *infile, *errorFile;
     int lastKUnitigNumber, kUnitigNumber=0, kUnitigNumberHold, i, argNum;
     int state;
     char *cptr, *cptr2, superReadName[32768], superReadNamePart[32768];
     char ori, oriHold;
     int overlap;
     int outputSeqLen;
     int seqDiffMax;
     int fail;
     char errorMessage[1000], errorMessageLine[1000];
     int numReads;
     char pluralStr[2];
     int noSequence = 0;
     
     // (VI) above
     workingDir = (char *) ".";
     sprintf (superReadListFile, "%s/%s", workingDir, DEFAULT_SUPER_READ_LIST_FILE);
     seqDiffMax = 0;
     errorFilename[0] = 0;
     argNum = 0;
     for (i=1; i<argc; i++) {
	  if (strcmp (argv[i], "-nosequence") == 0) {
	       noSequence = 1;
	       continue; }
	  if (strcmp (argv[i], "-seqdiffmax") == 0) {
	       ++i;
	       seqDiffMax = atoi (argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-error-filename") == 0) {
	       ++i;
	       strcpy (errorFilename, argv[i]);
	       continue; }
	  if (argNum == 0)
	       workingDir = argv[i];
	  else if (argNum == 1)
	       strcpy (superReadListFile, argv[i]);
	  else {
	       fprintf (stderr, "The program %s is called incorrectly. Too many args. Bye!\n", argv[0]);
	       exit (1); }
	  ++argNum; }
     
     if (strlen (errorFilename) == 0)
	  sprintf (errorFilename, "%s/createFastaSuperReadSequences.errors.txt", workingDir);
     mallocOrDie(line,MAX_READ_LEN, char);
     mallocOrDie (reverseComplementSpace, 1000000, char);
     mallocOrDie (outputSeqSpace, 3000000, char);

     sprintf (fname, "%s/%s", workingDir, KUNITIG_FILE_COMPLETE);
     stat (fname, &statbuf);
     fsize = statbuf.st_size;
     kUnitigSeqFileSize = fsize;
     // (I) above
     mallocOrDie (kUnitigSpace, fsize, char);
     // (II) above
     infile = Fopen (fname, "r");
     size_t bytes_read = fread (kUnitigSpace, 1, kUnitigSeqFileSize, infile);
     if(bytes_read != kUnitigSeqFileSize) {
       fprintf(stderr, "Failed to read the entire file '%s'. Bye!\n", fname);
       exit(2);
     }
     fclose (infile);

     // (III) above
     // Find out the last kUnitig number
     sprintf (kUnitigFilename, "%s/%s", workingDir, NUM_KUNITIGS_FILENAME);
     infile = Fopen (kUnitigFilename, "r");
     int fields_read = fscanf (infile, "%d\n", &lastKUnitigNumber);
     if(fields_read != 1) {
       fprintf(stderr, "Failed to read one int from '%s'. Bye!\n", kUnitigFilename);
       exit(2);
     }
     fclose (infile);
     fprintf (stderr, "The last kUnitigNumber was %d\n", lastKUnitigNumber);
     mallocOrDie (kUnitigSeq, lastKUnitigNumber+1, char *);
     mallocOrDie (kUnitigLengths, lastKUnitigNumber+1, int);

     // (IV) above
     state = AFTER_NEWLINE;
     for (i64=0; i64<kUnitigSeqFileSize; i64++) {
	  if (kUnitigSpace[i64] == '\n') {
	       kUnitigSpace[j64] = 0;
	       state = AFTER_NEWLINE;
	       continue; }
	  else if (state == IN_SEQUENCE) {
	       processTheChar; }
	  else if (state == AFTER_NEWLINE) {
	       if (kUnitigSpace[i64] == '>') {
		    ++i64;
		    kUnitigNumber = atoi (kUnitigSpace + i64);
		    while (kUnitigSpace[i64] != '\n')
			 ++i64;
		    state = AFTER_HEADER_NEWLINE;
		    continue; }
	       else { // It's after a regular newline
		    processTheChar; } }
	  else { // state == AFTER_HEADER_NEWLINE
	       j64 = i64;
	       kUnitigSeq[kUnitigNumber] = kUnitigSpace + i64;
	       processTheChar; }
     }

     for (i=0; i<=lastKUnitigNumber; i++)
	  if (kUnitigSeq[i] != NULL)
	       kUnitigLengths[i] = strlen (kUnitigSeq[i]);

#if 0
     for (i=0; i<=lastKUnitigNumber; i++)
	  if (kUnitigSeq[i] != NULL)
	       printf ("i = %d ; len = %d ; seq = %s\n", i, kUnitigLengths[i], kUnitigSeq[i]);
#endif

     // (V) above
     strcpy (fname, superReadListFile);
     infile = Fopen (fname, "r");
     errorFile = Fopen (errorFilename, "w");
     while (fgets (line, MAX_READ_LEN, infile)) {
	  numReads = 1;
	  cptr = line;
	  while ((cptr = strstr (cptr, " ; "))) {
	       ++numReads;
	       ++cptr; }
	  cptr = line;
	  while (isspace(*cptr)) ++cptr;
	  kUnitigNumber = atoi (cptr);
	  while (! isspace (*cptr)) ++cptr; while (isspace(*cptr)) ++cptr;
	  ori = *cptr;
	  kUnitigNumberHold = kUnitigNumber;
	  oriHold = ori;
	  sprintf (superReadName, "%d%c", kUnitigNumber, ori);
	  while (! isspace (*cptr)) ++cptr; while (isspace(*cptr)) ++cptr;
	  if (ori == 'F')
	       strcpy (outputSeqSpace, kUnitigSeq[kUnitigNumber]);
	  else {
	       generateReverseComplement (kUnitigSeq[kUnitigNumber], kUnitigLengths[kUnitigNumber]);
	       strcpy (outputSeqSpace, reverseComplementSpace); }
	  fail = 0;
	  errorMessage[0] = 0;
	  while (1) {
	       char *cptr3, *cptr4;
	       int numdiffs;
	       if (*cptr == ':')
		    break;
	       overlap = atoi (cptr);
	       while (! isspace (*cptr)) ++cptr; while (isspace(*cptr)) ++cptr;
	       kUnitigNumber = atoi (cptr);
	       while (! isspace (*cptr)) ++cptr; while (isspace(*cptr)) ++cptr;
	       ori = *cptr;
	       while (! isspace (*cptr)) ++cptr; while (isspace(*cptr)) ++cptr;
	       sprintf (superReadNamePart, "_%d_%d%c", overlap, kUnitigNumber, ori);
	       if (ori == 'F') {
		    cptr3 = kUnitigSeq[kUnitigNumber];
		    cptr2 = kUnitigSeq[kUnitigNumber] + overlap; }
	       else {
		    generateReverseComplement (kUnitigSeq[kUnitigNumber], kUnitigLengths[kUnitigNumber]);
		    cptr3 = reverseComplementSpace;
		    cptr2 = reverseComplementSpace + overlap; }
	       cptr4 = outputSeqSpace + (strlen(outputSeqSpace)-overlap);
	       numdiffs = 0;
	       for (i=0; i<overlap; i++)
		    if (cptr3[i] != cptr4[i])
			 ++numdiffs;
	       if (numdiffs > seqDiffMax) {
		    fail = 1;
		    if (numdiffs == 1)
			 strcpy (pluralStr, "");
		    else 
			 strcpy (pluralStr, "s");
		    sprintf (errorMessageLine, "The %d-base overlap between %d%c and %d%c has %d difference%s.\n", overlap, kUnitigNumberHold, oriHold, kUnitigNumber, ori, numdiffs, pluralStr);
		    strcat (errorMessage, errorMessageLine); }
	       strcat (outputSeqSpace, cptr2);
	       strcat (superReadName, superReadNamePart); 
	       kUnitigNumberHold = kUnitigNumber;
	       oriHold = ori; }
//     int outputSeqLen;
	  outputSeqLen = strlen (outputSeqSpace);
#if 0
	  printf ("superReadName = %s ; seq = %s\n", superReadName, outputSeqSpace);
#endif
	  if (fail) {
	       if (numReads == 1)
		    strcpy (pluralStr, "");
	       else
		    strcpy (pluralStr, "s");
	       sprintf (errorMessageLine, "%s (with %d read%s) fails\n", superReadName, numReads, pluralStr);
	       fputs (errorMessageLine, errorFile);
	       fputs (errorMessage, errorFile); fputc ('\n', errorFile);
	       continue; }
	  if (! noSequence)
	       fputc ('>', stdout);
	  fputs (superReadName, stdout); fputc ('\n', stdout);
          fflush(stdout);
	  if (! noSequence) {
	       fputs (outputSeqSpace, stdout); fputc ('\n', stdout); }
     }

     fclose (infile);
     fclose (errorFile);
	  
     return (0);
}

void generateReverseComplement (char *seq, int seqLen)
{
     int j, k;
     for (j=0, k=seqLen-1; k>=0; ++j, --k) {
	  switch (seq[j]) {
	  case 'a': reverseComplementSpace[k] = 't'; break;
	  case 'A': reverseComplementSpace[k] = 'T'; break;
	  case 'c': reverseComplementSpace[k] = 'g'; break;
	  case 'C': reverseComplementSpace[k] = 'G'; break;
	  case 'g': reverseComplementSpace[k] = 'c'; break;
	  case 'G': reverseComplementSpace[k] = 'C'; break;
	  case 't': reverseComplementSpace[k] = 'a'; break;
	  case 'T': reverseComplementSpace[k] = 'A'; break;
	  default: reverseComplementSpace[k] = seq[j]; break; }
     }
     reverseComplementSpace[seqLen] = 0;
}

FILE *Fopen (const char *fn, const char *mode)
{
     FILE *result;
     result = fopen (fn, mode);
     if (result == NULL)
     {
          fprintf (stderr, "Couldn't open file '%s' for ", fn);
          switch (mode[0]) {
          case 'r': fprintf (stderr, "reading"); break;
          case 'w': fprintf (stderr, "writing"); break;
          case 'a': fprintf (stderr, "appending"); break;
          default: fprintf (stderr, "unknown operation code '%c'", mode[0]);
               break;
          }
          fprintf (stderr, ". Bye!\n");
          exit (-1);
     }

     return (result);
}

