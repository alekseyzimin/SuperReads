#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <stdint.h>
#include <vector>
#include <charb.hpp>
#include <misc.hpp>

struct kUnitigOverlapStruct {
     int kUnitig1;
     int kUnitig2;
     int ahg;
     int bhg;
     char relativeOri;;
};

std::vector<char *> flds;
int *kUnitigLengths;
uint64_t *indexStarts;
int *kUnitigOvlVals;
unsigned char *howManyBasesInSpur;
unsigned char *kUniEndHasNonSpurFollow;
unsigned char *kUnitigIsDead;
unsigned char *kUnitigBeginOvlCount, *kUnitigEndOvlCount;

FILE *Fopen (const char *fn, const char *mode);

#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

int main (int argc, char **argv)
{
     char *workingDir = argv[1];
     int maxDeadEndBases = atoi (argv[2]);
//     int longestRgnToMerge = atoi (argv[3]);
     charb fname(1000);
     charb line(1000);
     FILE *infile;
     uint64_t numKUnitigs;
     int kmerLen = 31, kmerLenMinus1 = kmerLen-1;
     int tValue, maxBasesInDeadEndContinueForReadEnd;
     int numChanges, passNum, extraToPlayWith, tempSum;
     int numPossibleMerges=0;
     int killCount, numDeadForPass;
     kUnitigOverlapStruct kUnitigOverlapRecord;
     
     sprintf (fname, "%s/maxKUnitigNumber.txt", workingDir);
     infile = Fopen (fname, "r");
     fgets (line, 1000000, infile);
     numKUnitigs = atoll (line);
     fclose (infile);

     mallocOrDie (kUnitigLengths, numKUnitigs, int);
     sprintf (fname, "%s/kUnitigLengths.txt", workingDir);
     infile = Fopen (fname, "r");
     while (fgets (line, 1000000, infile)) {
	  getFldsFromLine (line, flds);
	  kUnitigLengths[atoll(flds[0])] = atoi (flds[1]); }
     fclose (infile);
     
     mallocOrDie (indexStarts, 2*numKUnitigs+1, uint64_t);
     
     sprintf (fname, "%s/overlap.overlaps", workingDir);
     infile = Fopen (fname, "rb");
     while (fread (&kUnitigOverlapRecord, sizeof (struct kUnitigOverlapStruct), 1, infile)) {
	  if (kUnitigOverlapRecord.ahg > 0) // Forward direction of read
	       ++indexStarts[2*kUnitigOverlapRecord.kUnitig1];
	  else
	       ++indexStarts[2*kUnitigOverlapRecord.kUnitig1+1];
     }
     for (uint64_t i=1; i<=2*numKUnitigs; i++)
	  indexStarts[i] += indexStarts[i-1];
     rewind (infile);
     mallocOrDie (kUnitigOvlVals, indexStarts[2*numKUnitigs], int);
     uint64_t tempIndex;
     while (fread (&kUnitigOverlapRecord, sizeof (struct kUnitigOverlapStruct), 1, infile)) {
	  if (kUnitigOverlapRecord.ahg > 0) // Forward direction of read
	       tempIndex = 2*kUnitigOverlapRecord.kUnitig1;
	  else
	       tempIndex = 2*kUnitigOverlapRecord.kUnitig1+1;
	  --indexStarts[tempIndex];
	  if (kUnitigOverlapRecord.relativeOri == 'N')
	       kUnitigOvlVals[indexStarts[tempIndex]] = 2*kUnitigOverlapRecord.kUnitig2;
	  else
	       kUnitigOvlVals[indexStarts[tempIndex]] = 2*kUnitigOverlapRecord.kUnitig2+1;
     }
     fclose (infile);

     mallocOrDie (howManyBasesInSpur, 2*numKUnitigs, unsigned char);
     numChanges = 1; // Needed some non-zero value to start
     passNum = 0;
     while (numChanges > 0) {
	  ++passNum;
//	  printf ("passNum = %d, numChanges = %d\n", passNum, numChanges);
	  numChanges = 0;
	  for (uint64_t i=0; i<2*numKUnitigs; i++) {
	       if (howManyBasesInSpur[i] != 0)
		    continue;
	       tValue = kUnitigLengths[i/2] - kmerLenMinus1;
	       if (tValue >= maxDeadEndBases)
		    continue;
	       extraToPlayWith = maxDeadEndBases - tValue;
	       if (indexStarts[i] == indexStarts[i+1]) {
		    if (tValue <= maxDeadEndBases) {
			 howManyBasesInSpur[i] = tValue;
			 ++numChanges; }
		    continue;
	       }
	       maxBasesInDeadEndContinueForReadEnd = 0;
	       for (uint64_t j=indexStarts[i]; j<indexStarts[i+1]; j++) {
		    int tempValue = kUnitigOvlVals[j], tempValue2;
		    tempSum = ((i%2) + (tempValue%2))%2;
		    tempValue2 = (tempValue/2)*2;
		    tempValue2 += tempSum;
		    if (howManyBasesInSpur[tempValue2] == 0) {
			 maxBasesInDeadEndContinueForReadEnd = maxDeadEndBases+1;
			 break; }
		    if (howManyBasesInSpur[tempValue2] >= extraToPlayWith) {
			 maxBasesInDeadEndContinueForReadEnd = maxDeadEndBases+1;
			 break; }
		    if (howManyBasesInSpur[tempValue2] >= maxBasesInDeadEndContinueForReadEnd)
			 maxBasesInDeadEndContinueForReadEnd = howManyBasesInSpur[tempValue2];
	       }
	       if (maxBasesInDeadEndContinueForReadEnd > maxDeadEndBases)
		    continue;
	       maxBasesInDeadEndContinueForReadEnd += tValue;
	       if (maxBasesInDeadEndContinueForReadEnd > maxDeadEndBases)
		    continue;
	       howManyBasesInSpur[i] = maxBasesInDeadEndContinueForReadEnd;
	       ++numChanges;
	  }
     }

//     printf ("howManyBasesInSpur:\n");
//     for (uint64_t i=0; i<numKUnitigs; i++) {
//	  printf ("%d %d %d\n", (int)i, howManyBasesInSpur[2*i], howManyBasesInSpur[2*i+1]); }
     mallocOrDie (kUniEndHasNonSpurFollow, 2*numKUnitigs, unsigned char);
     for (uint64_t i=0; i<2*numKUnitigs; i++) {
	  for (uint64_t j=indexStarts[i]; j<indexStarts[i+1]; j++) {
	       int tempValue = kUnitigOvlVals[j], tempValue2;
		    tempSum = ((i%2) + (tempValue%2))%2;
		    tempValue2 = (tempValue/2)*2;
		    tempValue2 += tempSum;
		    if (howManyBasesInSpur[tempValue2] == 0) {
			 kUniEndHasNonSpurFollow[i] = 1;
			 break; }
	  }
     }
//     printf ("howManyBasesInSpur:\n");
//     for (uint64_t i=0; i<numKUnitigs; i++) {
//	  printf ("kUniEndHasNonSpurFollow %d %d %d\n", (int)i, kUniEndHasNonSpurFollow[2*i], kUniEndHasNonSpurFollow[2*i+1]); }

     mallocOrDie (kUnitigIsDead, numKUnitigs, unsigned char);

     killCount = 0;
     for (uint64_t i=0; i<numKUnitigs; i++) {
	  int dontKill;
          //	  int path;
	  if ((howManyBasesInSpur[2*i] == 0) && (howManyBasesInSpur[2*i+1] == 0))
	       continue;
	  if (kUnitigIsDead[i])
	       continue;
	  dontKill = 0;
	  // Modify this so the spur to be cut can be shorter than what we need to
	  // consider the k-unitig extendable. Also change in following loop.
	  if (howManyBasesInSpur[2*i] > 0) { // Set for high-numbered bases
            //	       path = 1;
	       for (uint64_t j=indexStarts[2*i+1]; j<indexStarts[2*i+2]; j++) {
		    if (! kUniEndHasNonSpurFollow[kUnitigOvlVals[j]])
			 dontKill = 1;
		    if (dontKill)
			 break;
	       }
	       if (! dontKill) {
		    kUnitigIsDead[i] = 1;
		    if (kUnitigLengths[i] > 0)
			 ++killCount;
//		    if (i == 9804)
//			 printf ("At stage 1: k-uni = %d, path = %d\n", (int)i, path);
		    continue; } }
	  dontKill = 0;
	  if (howManyBasesInSpur[2*i+1] > 0) { // Set for high-numbered bases
            //	       path = 2;
	       for (uint64_t j=indexStarts[2*i]; j<indexStarts[2*i+1]; j++) {
		    int tempValue = kUnitigOvlVals[j];
		    tempValue = (tempValue/2)*2; // Make it even
		    if (kUnitigOvlVals[j] % 2 == 0)
			 ++tempValue;
		    if (! kUniEndHasNonSpurFollow[tempValue])
			 dontKill = 1;
		    if (dontKill)
			 break;
	       } 
	       if (! dontKill) {
		    kUnitigIsDead[i] = 1;
		    if (kUnitigLengths[i] > 0)
			 ++killCount;
//		    if (i == 9804)
//			 printf ("At stage 2: k-uni = %d, path = %d\n", (int)i, path);
		    continue; } }
     }
     printf ("%d k-unitigs killed\n", killCount);

     numDeadForPass = 1;
     while (numDeadForPass) {
	  int numOvlsNotDead, numOvlsDead;
	  numDeadForPass = 0;
//	  printf ("Starting a pass\n");
	  for (uint64_t i=0; i<numKUnitigs; i++) {
	       if (kUnitigIsDead[i])
		    continue;
	       if ((howManyBasesInSpur[2*i] == 0) && (howManyBasesInSpur[2*i+1] == 0))
		    continue;
	       numOvlsNotDead = numOvlsDead = 0;
	       for (uint64_t j=indexStarts[2*i]; j<indexStarts[2*i+1]; j++) {
		    if (kUnitigIsDead[kUnitigOvlVals[j]/2])
			 ++numOvlsDead;
		    else {
			 ++numOvlsNotDead;
			 break; } }
	       if ((numOvlsNotDead == 0) && (numOvlsDead != 0)) {
		    kUnitigIsDead[i] = 1;
//		    printf ("Killing k-unitig %d in part 1\n", (int)i);
		    ++numDeadForPass;
		    continue; }
	       if (kUnitigIsDead[i])
		    continue;
	       numOvlsNotDead = numOvlsDead = 0;
	       for (uint64_t j=indexStarts[2*i+1]; j<indexStarts[2*i+2]; j++) {
		    if (kUnitigIsDead[kUnitigOvlVals[j]/2])
			 ++numOvlsDead;
		    else {
			 ++numOvlsNotDead;
			 break; } }
	       if ((numOvlsNotDead == 0) && (numOvlsDead != 0)) {
		    kUnitigIsDead[i] = 1;
//		    printf ("Killing k-unitig %d in part 2\n", (int)i);
		    ++numDeadForPass;
		    continue; }
	  }
	  printf ("%d killed for the transitive dead-k-unitig killing\n", numDeadForPass);
	  killCount += numDeadForPass;
     }
     printf ("%d k-unitigs killed so far\n", killCount);
     printf ("Killed k-unitig list\n");
     for (uint64_t i=0; i<numKUnitigs; i++)
	  if (kUnitigIsDead[i] && kUnitigLengths[i])
	       printf ("%lld %d\n", (long long) i, kUnitigLengths[i]);
     mallocOrDie (kUnitigBeginOvlCount, numKUnitigs, unsigned char);
     mallocOrDie (kUnitigEndOvlCount, numKUnitigs, unsigned char);

     printf ("Num good joins for the k-unitigs that are left\n");
     for (uint64_t i=0; i<numKUnitigs; i++) {
	  if (kUnitigIsDead[i])
	       continue;
	  for (uint64_t j=indexStarts[2*i]; j<indexStarts[2*i+1]; j++)
	       if (! kUnitigIsDead[kUnitigOvlVals[j]/2])
		    ++kUnitigEndOvlCount[i];
	  for (uint64_t j=indexStarts[2*i+1]; j<indexStarts[2*i+2]; j++)
	       if (! kUnitigIsDead[kUnitigOvlVals[j]/2])
		    ++kUnitigBeginOvlCount[i];
	  printf ("%lld (len %d): %d %d\n", (long long) i, kUnitigLengths[i], kUnitigBeginOvlCount[i], kUnitigEndOvlCount[i]);
     }

     // Checking for merging
     for (uint64_t i=0; i<numKUnitigs; i++) {
	  if (kUnitigEndOvlCount[i] == 1) {
	       for (uint64_t j=indexStarts[2*i]; j<indexStarts[2*i+1]; j++) {
		    if (kUnitigIsDead[kUnitigOvlVals[j]/2])
			 continue;
		    // If we get here it's the unique overlap in this side
		    if (kUnitigOvlVals[j] % 2 == 0) {
			 if (kUnitigBeginOvlCount[kUnitigOvlVals[j]/2] == 1) {
			      ++numPossibleMerges; // Actually double counts
			      printf ("Merge %lld %d N 1\n", (long long) i, kUnitigOvlVals[j]/2); }
		    }
		    else {
			 if (kUnitigEndOvlCount[kUnitigOvlVals[j]/2] == 1) {
			      ++numPossibleMerges; // Actually double counts
			      printf ("Merge %lld %d I 1\n", (long long) i, kUnitigOvlVals[j]/2); }
		    }
		    break;
	       }
	  }
	  if (kUnitigBeginOvlCount[i] == 1) {
	       for (uint64_t j=indexStarts[2*i+1]; j<indexStarts[2*i+2]; j++) {
		    if (kUnitigIsDead[kUnitigOvlVals[j]/2])
			 continue;
		    // If we get here it's the unique overlap in this side
		    if (kUnitigOvlVals[j] % 2 == 0) {
			 if (kUnitigEndOvlCount[kUnitigOvlVals[j]/2] == 1) {
			      ++numPossibleMerges; // Actually double counts
			      printf ("Merge %lld %d N -1\n", (long long) i, kUnitigOvlVals[j]/2); }
		    }
		    else {
			 if (kUnitigBeginOvlCount[kUnitigOvlVals[j]/2] == 1) {
			      ++numPossibleMerges; // Actually double counts
			      printf ("Merge %lld %d I -1\n", (long long) i, kUnitigOvlVals[j]/2); }
		    }
		    break;
	       }
	  }
     }

     printf ("%d k-unitig pairs merged\n", numPossibleMerges/2);

     return (0);
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

