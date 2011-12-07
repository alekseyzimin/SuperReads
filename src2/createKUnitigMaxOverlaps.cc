/*
  The program outputs the coords and overlaps files for a given set of
  k-unitigs that overlap by exactly (k-1) bases. We work, by default,
  with k-unitigs generated using a k-mer size of 31, but if another
  value was used, use the flag
  -kmervalue kmerLen
  to specify the k-mer length used when generating the k-unitigs.

  The first non-flag arg is the prefix used for the k-unitigs file. (But see
  NOTE1 below.)
  The first flag assumes that the files are named *_#.fa, where
  * is the prefix specified and the #s start from 0 and continue until
  the last one.
  NOTE1: You may also use this arg to specify a filename. The filename must
  have the k-unitigs in numeric order.

  The second non-flag arg is the prefix used for the output files.
  The program will generate the files
  prefix.coords and prefix.overlaps.

  So the final syntax is
  createKUnitigMaxOverlaps [flags] inputPrefix outputPrefix OR
  createKUnitigMaxOverlaps [flags] inputFile   outputPrefix
  where the possible flags are
  -h: help and exit
  -kmervalue kmerLen
  -largest-kunitig-number largestKUnitigNumber

  NOTE2: This program assumes the k-unitig numbers are in increasing order
  in each file (unless the -largest-kunitig-number flag is used).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <assert.h>

#include <err.hpp>
#include <charb.hpp>

#define KMER_LENGTH 31
#define EST_OVLS_PER_KUNITIG 5

struct endKUnitigKmerStruct {
     unsigned long long kMerValue;
     int kUnitigNumber;
     unsigned char kUnitigEnd; // 0 or 1
     unsigned char ori; // 0 or 1
} *kMerMinusOneValuesAtEndOfKUnitigs, **ptrsToEndKUnitigKmerStructs;

struct overlapDataStruct
{
     int kUni1;
     int kUni2; 
     int ahg;
     int bhg;
     char netOri;
} *overlapData;

charb line(2000);
char **kUnitigSequences;
unsigned char **kUnitigSequenceCounts;
unsigned char *endIsDone;
int *kUnitigLengths, largestKUnitigNumber;
int kmerLen;
char *inputPrefix, *outputPrefix;
int *startOverlapByUnitig;
struct overlapDataStruct *overlapDataToSave;

void reportKUnitigEndMatches (void);
int kmerStructCompare (struct endKUnitigKmerStruct **ptr1, struct endKUnitigKmerStruct **ptr2);
void loadKUnitigEndingKMerValues (void);
int getLargestKUnitigNumber (char *prefix, int numInputFiles);
int getNumInputFiles (char *prefix);
void loadKUnitigSequences (char *prefix, int numInputFiles);
FILE *Fopen (const char *fn, const char *mode);
void giveUsageAndExit (void);
void processArgs (int argc, char **argv);

#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

int main (int argc, char *argv[])
{
     int numInputFiles, i;

     processArgs (argc, argv);
     numInputFiles = getNumInputFiles (inputPrefix);
     if (numInputFiles >= 1)
	  printf ("num input files = %d\n", numInputFiles);
     else
	  printf ("input file = %s\n", inputPrefix);

     if (largestKUnitigNumber == 0)
	  largestKUnitigNumber = getLargestKUnitigNumber (inputPrefix, numInputFiles);
     printf ("largestKUnitigNumber = %d\n", largestKUnitigNumber);

     mallocOrDie (kUnitigSequences, largestKUnitigNumber+1, char *);
     mallocOrDie (kUnitigLengths, largestKUnitigNumber+1, int);
     mallocOrDie (overlapData, (largestKUnitigNumber+1) * EST_OVLS_PER_KUNITIG, struct overlapDataStruct);
     mallocOrDie (startOverlapByUnitig, largestKUnitigNumber+2, int);
     loadKUnitigSequences (inputPrefix, numInputFiles);

     mallocOrDie (kMerMinusOneValuesAtEndOfKUnitigs, 4*(largestKUnitigNumber+1), struct endKUnitigKmerStruct);
     loadKUnitigEndingKMerValues ();
     mallocOrDie (ptrsToEndKUnitigKmerStructs, 4*(largestKUnitigNumber+1), struct endKUnitigKmerStruct *);
     for (i=0; i<4*(largestKUnitigNumber+1); i++)
	  ptrsToEndKUnitigKmerStructs[i] = kMerMinusOneValuesAtEndOfKUnitigs+i;
     qsort (ptrsToEndKUnitigKmerStructs, 4*(largestKUnitigNumber+1), sizeof (struct endKUnitigKmerStruct *), (int (*)(const void*, const void*)) kmerStructCompare);

     mallocOrDie (endIsDone, largestKUnitigNumber+1, unsigned char);

     reportKUnitigEndMatches ();

     return (0);
}

void reportKUnitigEndMatches (void)
{
     int beginIndex, endIndex;
     struct endKUnitigKmerStruct *ptr1, *ptr2;
     int i, j, totKUniStartSep;
     int kUni1, kUni2, ahg, bhg, begin1, end1, begin2, end2;
     int isGoodOverlap, skipThis;
     uint64_t numOvlsOutput=0;
     char netOri;
     char filename[500];
     FILE *coordsFile, *overlapsFile;

     sprintf (filename, "%s.coords", outputPrefix);
     coordsFile = Fopen (filename, "w");
     sprintf (filename, "%s.overlaps", outputPrefix);
     overlapsFile = Fopen (filename, "wb");

     beginIndex=0;
     while (beginIndex<4*(largestKUnitigNumber+1)) {
	  ptr1 = ptrsToEndKUnitigKmerStructs[beginIndex];
	  for (endIndex=beginIndex+1; (endIndex<4*(largestKUnitigNumber+1)) && (ptrsToEndKUnitigKmerStructs[endIndex]->kMerValue == ptr1->kMerValue); endIndex++) {
	  }
	  skipThis = 0;
	  if (ptrsToEndKUnitigKmerStructs[beginIndex]->kUnitigEnd == 0) {
	       if (endIsDone[ptrsToEndKUnitigKmerStructs[beginIndex]->kUnitigNumber] >= 2)
		    skipThis = 1;
	  }
	  else {
	       if (endIsDone[ptrsToEndKUnitigKmerStructs[beginIndex]->kUnitigNumber] % 2 == 1)
		    skipThis = 1;
	  }
	  if (skipThis)
	       goto endOfLoop;
	  for (i=beginIndex; i<endIndex; i++) {
	       ptr1 = ptrsToEndKUnitigKmerStructs[i];
	       kUni1 = ptr1->kUnitigNumber;
	       if (ptr1->kUnitigEnd == 0)
		    endIsDone[kUni1] += 2;
	       else
		    endIsDone[kUni1] += 1;
	       if (kUnitigLengths[kUni1] == 0)
		    continue;
	       
	       for (j=beginIndex; j<endIndex; j++) {
		    if (i == j)
			 continue;
		    ptr2 = ptrsToEndKUnitigKmerStructs[j];
		    kUni2 = ptr2->kUnitigNumber;
		    if (kUnitigLengths[kUni2] == 0)
			 continue;
		    if (ptr1->ori == ptr2->ori)
			 netOri = 'N';
		    else
			 netOri = 'I';
		    if ((ptr1->kUnitigEnd + ptr2->kUnitigEnd + ptr1->ori + ptr2->ori) % 2 == 0)
			 isGoodOverlap = 0;
		    else
			 isGoodOverlap = 1;
		    if (isGoodOverlap) {
			 if (ptr1->kUnitigEnd == 0) {
			      ahg = (kmerLen-1) - kUnitigLengths[kUni2];
			      bhg = (kmerLen-1) - kUnitigLengths[kUni1]; }
			 else {
			      ahg = kUnitigLengths[kUni1] - (kmerLen-1);
			      bhg = kUnitigLengths[kUni2] - (kmerLen-1); }
		    }
		    else {
			 if (ptr1->kUnitigEnd == 0) {
			      ahg = 0;
			      bhg = kUnitigLengths[kUni2] - kUnitigLengths[kUni1]; }
			 else {
			      ahg = kUnitigLengths[kUni1] - kUnitigLengths[kUni2];
			      bhg = 0; }
		    }
		    if (ptr1->kUnitigEnd == 0) {
			 begin1 = 1;
			 end1 = kmerLen-1; }
		    else {
			 end1 = kUnitigLengths[kUni1];
			 begin1 = end1 - (kmerLen-1) + 1; }
		    if (netOri == 'N') {
			 begin2 = begin1 - ahg;
			 end2 = end1 - ahg; }
		    else {
			 totKUniStartSep = kUnitigLengths[kUni1] + bhg;
			 begin2 = totKUniStartSep - (begin1-1);
			 end2 = (totKUniStartSep - end1) + 1;
		    }
		    fprintf (coordsFile, "%d %d %d %d 100.00 %d %d %d %d\n", begin1, end1, begin2, end2, kUnitigLengths[kUni1], kUnitigLengths[kUni2], kUni1, kUni2);
		    if (isGoodOverlap) {
			 if(kUni1 !=kUni2){//temporary dirty fix by Aleksey
			      overlapData[numOvlsOutput].kUni1 = kUni1;
			      overlapData[numOvlsOutput].kUni2 = kUni2;
			      overlapData[numOvlsOutput].ahg = ahg;
			      overlapData[numOvlsOutput].bhg = bhg;
			      overlapData[numOvlsOutput].netOri = netOri;
			      ++numOvlsOutput;
//			 fprintf (overlapsFile, "%d %d %c %d %d 0.0 0.0\n", kUni1, kUni2, netOri, ahg, bhg);
			 }
		    }
	       }
	  }
     endOfLoop:
	  beginIndex = endIndex;
     }

     mallocOrDie (overlapDataToSave, numOvlsOutput, struct overlapDataStruct);
     for (uint64_t j=0; j<numOvlsOutput; j++)
     {
	  int unitig1 = overlapData[j].kUni1, unitig2 = overlapData[j].kUni2;
	  
#if 0
	  if (unitig1 > unitig2) 
	       continue;
	  else if ((unitig1 == unitig2) && (overlapData[j].ahg < 0))
	       continue;
#endif
	  if (unitig1 >= unitig2)
	       continue;
	  startOverlapByUnitig[unitig1]++;
	  startOverlapByUnitig[unitig2]++;
     }
     for (int64_t unitigNum = 1; unitigNum < largestKUnitigNumber + 2; unitigNum++)
	  startOverlapByUnitig[unitigNum] += startOverlapByUnitig[unitigNum - 1];
     
     for (uint64_t j=0; j<numOvlsOutput; j++)
     {
	  int unitig1 = overlapData[j].kUni1, unitig2 = overlapData[j].kUni2;
	  int ahg = overlapData[j].ahg, bhg = overlapData[j].bhg;
	  char ori = overlapData[j].netOri;
	  int itemp, itempHold;
	  
#if 0
	  if (unitig1 > unitig2) 
	       continue;
	  else if ((unitig1 == unitig2) && (ahg < 0))
	       continue;
#endif
	  if (unitig1 >= unitig2)
	       continue;
	  startOverlapByUnitig[unitig1]--;
	  itemp = startOverlapByUnitig[unitig1];
	  overlapDataToSave[itemp].kUni1 = unitig1;
	  overlapDataToSave[itemp].kUni2 = unitig2;
	  overlapDataToSave[itemp].netOri = ori;
	  overlapDataToSave[itemp].ahg = ahg;
	  overlapDataToSave[itemp].bhg = bhg;
	  startOverlapByUnitig[unitig2]--;
	  itempHold = itemp;
	  itemp = startOverlapByUnitig[unitig2];
	  overlapDataToSave[itemp].kUni1 = unitig2;
	  overlapDataToSave[itemp].kUni2 = unitig1;
	  overlapDataToSave[itemp].netOri = ori;
	  if (ori == 'N')
	  {
	       overlapDataToSave[itemp].ahg = -ahg;
	       overlapDataToSave[itemp].bhg = -bhg;
	  }
	  else
	  {
	       overlapDataToSave[itemp].ahg = bhg;
	       overlapDataToSave[itemp].bhg = ahg;
	  }
     }

     fwrite (overlapDataToSave, sizeof (struct overlapDataStruct), numOvlsOutput, overlapsFile);
     
     fclose (coordsFile);
     fclose (overlapsFile);
}
	  
int kmerStructCompare (struct endKUnitigKmerStruct **ptr1, struct endKUnitigKmerStruct **ptr2)
{
     if ((*ptr1)->kMerValue < (*ptr2)->kMerValue) return (-1);
     if ((*ptr1)->kMerValue > (*ptr2)->kMerValue) return (1);
     if ((*ptr1)->kUnitigNumber < (*ptr2)->kUnitigNumber) return (-1);
     if ((*ptr1)->kUnitigNumber > (*ptr2)->kUnitigNumber) return (1);
     if ((*ptr1)->kUnitigEnd < (*ptr2)->kUnitigEnd) return (-1);
     if ((*ptr1)->kUnitigEnd > (*ptr2)->kUnitigEnd) return (1);
     if ((*ptr1)->ori < (*ptr2)->ori) return (-1);
     if ((*ptr1)->ori > (*ptr2)->ori) return (1);
     return (0);
}

void loadKUnitigEndingKMerValues (void)
{
     unsigned long long mask = 0ULL;
     int i, kUnitigNumber;
     unsigned long index;
     
     for (i=0; i<32; i++) {
	  mask <<= 2;
	  mask += 3; }
     for (kUnitigNumber=0; kUnitigNumber<=largestKUnitigNumber; kUnitigNumber++) {
	  index = 4 * kUnitigNumber;
	  if (kUnitigLengths[kUnitigNumber] == 0) {
	       kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue = kMerMinusOneValuesAtEndOfKUnitigs[index+1].kMerValue = kMerMinusOneValuesAtEndOfKUnitigs[index+2].kMerValue = kMerMinusOneValuesAtEndOfKUnitigs[index+3].kMerValue = mask;
	       kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigNumber = kMerMinusOneValuesAtEndOfKUnitigs[index+1].kUnitigNumber = kMerMinusOneValuesAtEndOfKUnitigs[index+2].kUnitigNumber = kMerMinusOneValuesAtEndOfKUnitigs[index+3].kUnitigNumber = kUnitigNumber;
	       kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigEnd = kMerMinusOneValuesAtEndOfKUnitigs[index+1].kUnitigEnd = 0;
	       kMerMinusOneValuesAtEndOfKUnitigs[index+2].kUnitigEnd = kMerMinusOneValuesAtEndOfKUnitigs[index+3].kUnitigEnd = 1;
	       kMerMinusOneValuesAtEndOfKUnitigs[index].ori = kMerMinusOneValuesAtEndOfKUnitigs[index+2].ori = 0;
	       kMerMinusOneValuesAtEndOfKUnitigs[index+1].ori = kMerMinusOneValuesAtEndOfKUnitigs[index+3].ori = 1;
	       continue;
	  }
	  // k-mer at beginning of k-unitig, forward ori
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue = 0ULL;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigNumber = kUnitigNumber;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigEnd = 0;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].ori = 0;
	  for (i=0; i<kmerLen-1; i++) {
	       kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue <<= 2;
	       switch (kUnitigSequences[kUnitigNumber][i]) {
	       case 'a': case 'A': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 0; break;
	       case 'c': case 'C': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 1; break;
	       case 'g': case 'G': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 2; break;
	       case 't': case 'T': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 3; break;
	       }
	  }
	  ++index;
	  // k-mer at beginning of k-unitig, reverse ori
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue = 0ULL;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigNumber = kUnitigNumber;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigEnd = 0;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].ori = 1;
	  for (i=kmerLen-2; i>=0; i--) {
	       kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue <<= 2;
	       switch (kUnitigSequences[kUnitigNumber][i]) {
	       case 't': case 'T': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 0; break;
	       case 'g': case 'G': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 1; break;
	       case 'c': case 'C': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 2; break;
	       case 'a': case 'A': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 3; break;
	       }
	  }
	  ++index;
	  // k-mer at end of k-unitig, forward ori
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue = 0ULL;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigNumber = kUnitigNumber;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigEnd = 1;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].ori = 0;
	  for (i=kUnitigLengths[kUnitigNumber]-kmerLen+1; i<kUnitigLengths[kUnitigNumber]; i++) {
	       kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue <<= 2;
	       switch (kUnitigSequences[kUnitigNumber][i]) {
	       case 'a': case 'A': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 0; break;
	       case 'c': case 'C': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 1; break;
	       case 'g': case 'G': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 2; break;
	       case 't': case 'T': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 3; break;
	       }
	  }
	  ++index;
	  // k-mer at end of k-unitig, reverse ori
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue = 0ULL;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigNumber = kUnitigNumber;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].kUnitigEnd = 1;
	  kMerMinusOneValuesAtEndOfKUnitigs[index].ori = 1;
	  for (i=kUnitigLengths[kUnitigNumber]-1; i>kUnitigLengths[kUnitigNumber]-kmerLen; i--) {
	       kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue <<= 2;
	       switch (kUnitigSequences[kUnitigNumber][i]) {
	       case 't': case 'T': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 0; break;
	       case 'g': case 'G': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 1; break;
	       case 'c': case 'C': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 2; break;
	       case 'a': case 'A': kMerMinusOneValuesAtEndOfKUnitigs[index].kMerValue += 3; break;
	       }
	  }
     }
}	       

int getLargestKUnitigNumber (char *prefix, int numInputFiles)
{
     struct stat statbuf;
     FILE *infile;
     int i;
     char fname[500];
     off_t fileOffset;
     unsigned int maxKUnitig = 0;
     unsigned int currentKUnitig;
     int numInputFilesTemp;

     numInputFilesTemp = (numInputFiles == 0) ? 1 : numInputFiles;

     for (i=0; i<numInputFilesTemp; i++) {
	  if (numInputFiles > 0)
	       sprintf (fname, "%s_%d.fa", prefix, i);
	  else
	       strcpy (fname, prefix);
	  stat (fname, &statbuf);
	  fileOffset = statbuf.st_size;
	  if (fileOffset >= 1000000)
	       fileOffset -= 1000000;
	  else
	       fileOffset = 0;
	  infile = Fopen (fname, "r");
	  fseek (infile, fileOffset, SEEK_SET);
	  if(!fgets (line, infile))
            die << "Error reading file '" << fname << "'" << err::no;
	  while (fgets (line, 1000000, infile)) {
	       if (line[0] != '>')
		    continue;
	       currentKUnitig = atoi (line+1);
	       if (currentKUnitig > maxKUnitig)
		    maxKUnitig = currentKUnitig;
	  }
	  fclose (infile);
     }
     return (maxKUnitig);
}

int getNumInputFiles (char *prefix)
{
     struct stat statbuf;
     int i;
     char fname[500];
     for (i=0; 1; i++) {
	  sprintf (fname, "%s_%d.fa", prefix, i);
	  if (stat (fname, &statbuf) == 0)
	       continue;
	  return (i);
     }
     return (0);
}

void loadKUnitigSequences (char *prefix, int numInputFiles)
{
     FILE *infile;
     int i, kUnitigNumber, length;
     char fname[500];
     int numInputFilesTemp;
     charb header;

     numInputFilesTemp = (numInputFiles == 0) ? 1 : numInputFiles;

     for (i=0; i<numInputFilesTemp; i++) {
	  if (numInputFiles > 0)
	       sprintf (fname, "%s_%d.fa", prefix, i);
	  else
	       strcpy (fname, prefix);
	  infile = Fopen (fname, "r");
	  if(!infile)
	    die << "Can't open file '" << fname << "'" << err::no;
	  
	  int next_char = fgetc(infile);
	  if(next_char != '>')
	    die << "Badly formatted fasta file '" << fname << "'. Missing header";

	  while(fgets (header, infile)) {
	    kUnitigNumber = atoi (header);
	    next_char = fgetc(infile);
	    line.clear();
	    while(next_char != EOF && next_char != '>') {
	      ungetc(next_char, infile);
	      fgets_append (line, infile);
	      fflush(stdout);
	      line.chomp();
	      next_char = fgetc(infile);
	    }
	    length = line.len();
	    mallocOrDie (kUnitigSequences[kUnitigNumber], length+1, char);
	    kUnitigLengths[kUnitigNumber] = length;
	    strcpy (kUnitigSequences[kUnitigNumber], line);
	  }
	  fclose (infile);
     }
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

void giveUsageAndExit (void)
{
     fprintf (stderr, "This program outputs the coords and overlaps files for a given set of\n k-unitigs generated by k-mers of length K that overlap by exactly\n (K-1) bases. We work, by default, with k-unitigs generated using a\n k-mer size of 31, but if another k-mer size was used, use the flag\n-kmervalue kMerSize\n to specify the k-mer size used when generating the k-unitigs.\n\nThe first non-flag arg is the prefix used for the k-unitigs files.\n It assumes that the files are named *_#.fa, where\n * is the prefix specified and the #s start from 0 and continue until\n the last one. This arg may also be used to specify the complete\n filename. Note that all input files are assumed to have k-unitig\n numbers in ascending order.\n\nThe second non-flag arg is the prefix used for the output files.\n The program will generate the files\n prefix.coords   and   prefix.overlaps.\n\n So the final syntax is\n\ncreateKUnitigMaxOverlaps [flags] inputPrefix outputPrefix\n where the possible flags are\n   -h: help and exit\n   -kmervalue kMerSize\n   -largest-kunitig-number largestKUnitigNumber (in this case the\n       k-unitigs don't have to be in numeric order in the files.)\n");
     exit (0);
}

void processArgs (int argc, char **argv)
{
     int numArgsSeen, i;
     inputPrefix = outputPrefix = NULL;
     kmerLen = KMER_LENGTH;
     numArgsSeen = 0;
     largestKUnitigNumber = 0;
     for (i=1; i<argc; i++) {
	  if (strcmp (argv[i], "-h") == 0)
	       giveUsageAndExit();
	  if (strcmp (argv[i], "-kmervalue") == 0) {
	       ++i;
	       kmerLen = atoi (argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-largest-kunitig-number") == 0) {
	       ++i;
	       largestKUnitigNumber = atoi (argv[i]);
	       continue; }
	  if (argv[i][0] == '-') {
	       fprintf (stderr, "\nUnrecognized flag: %s\n\n", argv[i]);
	       giveUsageAndExit(); }
	  if (numArgsSeen == 0)
	       inputPrefix = argv[i];
	  else if (numArgsSeen == 1)
	       outputPrefix = argv[i];
	  else {
	       fprintf (stderr, "\nToo many args supplied. 2 are expected.\n\n");
	       giveUsageAndExit(); }
	  ++numArgsSeen;
     }
     if (numArgsSeen < 2) {
	  fprintf (stderr, "\nToo few args supplied. 2 are expected.\n\n");
	  giveUsageAndExit();
     }
}

