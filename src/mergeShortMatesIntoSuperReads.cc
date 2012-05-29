/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <libgen.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <src/diskBasedUnitigger.h>
#include "charb.hpp"
using namespace std;

#define OFFSET_MEANING_END 126
#define OFFSET_MEANING_END_DOUBLED (OFFSET_MEANING_END << 1)
#define MAX_ARRAY_ELEMENTS_TEMP 1000000000

struct superReadStruct {
     unsigned int *kUnitigList;
     unsigned char *offsetAndOriList; // offset<<1+ori; if (offset ==126) is end
} srsData1, srsData2, srsData3, srsData4, badSrsData, *badSrsDataSentinel, srsGlobalData, *srsResultsPtr;

int minSepNonRptKUni;
char superReadsFile[1000], numKUnitigsFile[1000], chimericReadsFile[1000], repetitiveKUnitigsFile[1000], kUnitigLengthsFile[1000], statsFile[1000], numReadsFile[1000];
int *kUnitigLengths;
unsigned int *kUnitigList1;

char *flds[1000];
char forwardReadLine2[1000], reverseReadLine2[1000];
unsigned char *isRepetitive, *isChimeric, *wasOutput, *insertCount, *insertEndForFirstReadOfInsert;
unsigned int numKUnitigs;
unsigned long long numReads, numInserts;
unsigned long long *startIndex;
unsigned long long identicalMateSuperReads, numContainedSuperReadsAt1, numContainedSuperReadsAt2, newlyMerged, forcedMerges, diffAndFailToMerge;

void processArgs (int argc, char *argv[]);
void setStatsFilename (char *thisProgramName);
int getFldsFromLine (char *cptrHold);
int getSuperReadLengthFromSuperReadStruct (struct superReadStruct *srsPtr);
int getNumKUnitigsInSuperReadStruct (struct superReadStruct *srsPtr);
// In the following we assume the space for the struct outputStsPtr is already
// allocated
void returnReversedSuperReadData (struct superReadStruct *inputSrsPtr, struct superReadStruct *outputSrsPtr);
struct superReadStruct *getMinOfSuperReadAndReverseComplement (struct superReadStruct *inputSrsPtr);
// kUnitigList in next function expects a sentinel at end >= numKUnitigs
int reportIfSomeKUnitigInArrayIsNonRepetitive (unsigned int *kUnitigList);
// The following assumes the appropriate memory has been allocated.
// The kUnitigList ends in a sentinel >= numKUnitigs
void getKUnitigsFromSuperRead (struct superReadStruct *srsPtr, unsigned int *kUnitigList);
// The return value is how many k-unitigs are in the super-read
int setSrsStruct (struct superReadStruct *srsInput, unsigned long long indexIn, struct superReadStruct *srsOutput, unsigned long long indexOut);
void printSrsStruct (char *str, struct superReadStruct *srsPtr);
struct superReadStruct *calculateInsertSuperRead (struct superReadStruct *inputSrsPtr1, struct superReadStruct *inputSrsPtr2, struct superReadStruct *outputSrsPtr);
int superReadIsReversedForNormalForm (struct superReadStruct *srsPtr);
void returnSuperReadInNormalForm (struct superReadStruct *inputSrsPtr, struct superReadStruct *outputSrsPtr);

#define getOffsetAndOriFromCode(code, overlap, ori) ori = code % 2; overlap = code >> 1
#define shiftCodeAndIfAtEnd(code) code >>= 1; if (code >= OFFSET_MEANING_END)
#define outputReadLines(readName, readSrsData, outfile) fprintf (outfile, "readName = %llu : ", readName); \
printSrsStruct (outputLine, &readSrsData); \
fputs (outputLine, outfile); \
fputc ('\n', outfile); \
if (superReadIsReversedForNormalForm (&readSrsData)) { \
returnReversedSuperReadData (&readSrsData, &srsData3); \
printSrsStruct (outputLine, &srsData3); } \
fputs (outputLine, outfile); \
fputc ('\n', outfile)
#define initializeSrsData(srsName, numElements) mallocOrDie (srsName.kUnitigList, numElements, unsigned int); \
     mallocOrDie (srsName.offsetAndOriList, numElements, unsigned char);
#define READ_NUM_TO_STUDY 1021713
// #define DE if (readName == READ_NUM_TO_STUDY)
#define DE if (0)
#define OUTPUT_FILE stdout

int main (int argc, char **argv)
{
     char line2[1000000], outputLine[1000000];
     charb line(10);
     unsigned char tOri, tOri2, ucChar2;
     int numFlds;
     long long lltemp;
     unsigned long long currentIndex=0LLU, ullTemp;
     FILE *infile, *outfile;
     int itemp;
     int j, k;
     struct superReadStruct fwdSrsData, revSrsData, reversedRevSrsData;
     struct superReadStruct reversedFwdSrsData;
     int someKUnitigIsNotRepetitive, isSubstring;
     int numFwdFields, numRevFields;
     int forwardLen, reverseLen;

     processArgs (argc, argv);
     badSrsDataSentinel = &badSrsData;
     identicalMateSuperReads = numContainedSuperReadsAt1 = numContainedSuperReadsAt2 = newlyMerged = forcedMerges = diffAndFailToMerge = 0LLU;
     setStatsFilename (argv[0]);
//     mallocOrDie (srsData1.kUnitigList, 10, unsigned int);
//     mallocOrDie (srsData1.offsetAndOriList, 10, unsigned char);
     initializeSrsData (srsData1, 1000);
     initializeSrsData (srsData2, 1000);
     initializeSrsData (srsData3, 1000);
     initializeSrsData (srsData4, 1000);
     initializeSrsData (fwdSrsData, 1000);
     initializeSrsData (revSrsData, 1000);
     initializeSrsData (reversedRevSrsData, 1000);
     initializeSrsData (reversedFwdSrsData, 1000);
     mallocOrDie (kUnitigList1, 1000, unsigned int);
     mallocOrDie (startIndex, MAX_ARRAY_ELEMENTS_TEMP, unsigned long long);
     mallocOrDie (srsGlobalData.kUnitigList, 2*MAX_ARRAY_ELEMENTS_TEMP, unsigned int);
     mallocOrDie (srsGlobalData.offsetAndOriList, 2*MAX_ARRAY_ELEMENTS_TEMP, unsigned char);

     fprintf (stderr, "Opening the numKUnitigsFile...\n");
     infile = Fopen (numKUnitigsFile, "r");
     fscanf (infile, "%u\n", &numKUnitigs);
     fclose (infile);

     mallocOrDie (isRepetitive, numKUnitigs+1, unsigned char);

     fprintf (stderr, "Opening the repetitiveKUnitigsFile...\n");
     infile = Fopen (repetitiveKUnitigsFile, "r");
     while (fgets (line, 1000000, infile)) {
	  numFlds = getFldsFromLine (line);
	  itemp = atoi (flds[1]);
	  if (itemp >= minSepNonRptKUni)
	       continue;
	  isRepetitive[atoi(flds[0])] = 1;
     }
     fclose (infile);

     fprintf (stderr, "Opening the numReadsFile...\n");
     infile = Fopen (numReadsFile, "r");
     fscanf (infile, "%llu\n", &numReads);
     fclose (infile);

     numInserts = numReads/2;

     mallocOrDie (isChimeric, numReads+1, unsigned char);

     fprintf (stderr, "Opening the chimericReadsFile...\n");
     infile = Fopen (chimericReadsFile, "r");
     while (fgets (line, 1000000, infile)) {
	  lltemp = atoll (line);
	  isChimeric[lltemp] = 1;
     }
     fclose (infile);
     
     mallocOrDie (kUnitigLengths, numKUnitigs+1, int);
     infile = Fopen (kUnitigLengthsFile, "r");
     while (fgets (line, 1000000, infile)) {
	  numFlds = getFldsFromLine (line);
	  kUnitigLengths[atoi(flds[0])] = atoi (flds[1]);
     }
     fclose (infile);

     mallocOrDie (wasOutput, numReads+1, unsigned char);
     mallocOrDie (insertCount, numInserts+1, unsigned char);
     mallocOrDie (insertEndForFirstReadOfInsert, numInserts+1, unsigned char);
     fprintf (stderr, "Opening the superReadsFile...\n");
     infile = Fopen (superReadsFile, "r");
     while (fgets (line, 1000000, infile)) {
	  char *cptr;
	  unsigned long long readName, otherRead, readInsertNum;
	  unsigned long long fwdReadName, revReadName;
	  unsigned char readInsertEnd;
//	  fprintf (stderr, "We got to 1\n");
	  fgets (line2, 1000000, infile);
	  cptr = strstr ((char *)line, "readName = ");
	  cptr += strlen ("readName = ");
	  sscanf (cptr, "%llu", &readName);
	  DE fprintf (OUTPUT_FILE, "We got to 2\n");
	  if (isChimeric[readName]) {
	       wasOutput[readName] = 1;
	       continue; }
	  DE fprintf (OUTPUT_FILE, "We got to 3\n");
	  cptr = strchr (cptr, ':');
	  cptr += 2;
	  readInsertNum = readName / 2;
	  readInsertEnd = readName % 2;
	  DE fprintf (OUTPUT_FILE, "We got to 10\n");

	  numFlds = getFldsFromLine (cptr);
	  for (j=0; j<numFlds; j+=3) {
	       int index = j/3;
	       srsData1.kUnitigList[index] = atoi (flds[j]);
	       if (j+2<numFlds)
		    srsData1.offsetAndOriList[index] = atoi (flds[j+2]);
	       else
		    srsData1.offsetAndOriList[index] = OFFSET_MEANING_END;
	       (srsData1.offsetAndOriList)[index] <<= 1;
	       if (flds[j+1][0] == 'R')
		    (srsData1.offsetAndOriList)[index]+=1;
	       DE printf ("index = %d, unitig = %d, offsetAndOri = %d\n", index, srsData1.kUnitigList[index], (srsData1.offsetAndOriList)[index]);
	  }
	  DE fprintf (OUTPUT_FILE, "We got to 20\n");
	  k = (j+1)/3;
	  ++insertCount[readInsertNum];
	  DE fprintf (OUTPUT_FILE, "readInsertNum = %llu, insertCount = %d\n", readInsertNum, insertCount[readInsertNum]);
	  if (insertCount[readInsertNum] == 1) {
	       insertEndForFirstReadOfInsert[readInsertNum] = readInsertEnd;
	       startIndex[readInsertNum] = currentIndex; // Nothing defined or initialized for this yet
	       k = setSrsStruct (&srsData1, 0, &srsGlobalData, currentIndex);
	       currentIndex += k;
	       continue;
	  }
	  DE fprintf (OUTPUT_FILE, "We got to 30\n");
	  // If we get here it's the 2nd read for the insert
	  if (readName % 2 == 0) {
	       otherRead = readName + 1;
	       fwdReadName = readName;
	       numFwdFields = setSrsStruct (&srsData1, 0, &fwdSrsData, 0);
	       numRevFields = setSrsStruct (&srsGlobalData, startIndex[readInsertNum], &revSrsData, 0);
	  }
	  else {
	       otherRead = readName - 1;
	       fwdReadName = otherRead;
	       numFwdFields = setSrsStruct (&srsGlobalData, startIndex[readInsertNum], &fwdSrsData, 0);
	       numRevFields = setSrsStruct (&srsData1, 0, &revSrsData, 0);
	  }
	  DE printf ("otherRead = %llu, readInsertNum = %llu, startIndex[readInsertNum] = %llu\n", otherRead, readInsertNum, startIndex[readInsertNum]);
	  DE fprintf (OUTPUT_FILE, "We got to 40\n");
	  revReadName = fwdReadName + 1;
	  returnReversedSuperReadData (&revSrsData, &reversedRevSrsData);
	  DE fprintf (OUTPUT_FILE, "We got to 45, numFwdFields = %d, numRevFields = %d\n", numFwdFields, numRevFields);
	  if (numFwdFields == numRevFields) {
	       int isEqual = 1;
	       for (k=0; k<numFwdFields; k++) {
		    if (fwdSrsData.kUnitigList[k] != reversedRevSrsData.kUnitigList[k])
			 isEqual = 0;
		    if (fwdSrsData.offsetAndOriList[k] != reversedRevSrsData.offsetAndOriList[k])
			 isEqual = 0;
		    if (! isEqual)
			 break;
	       }
	       DE fprintf (OUTPUT_FILE, "We got to 48\n");
	       if (isEqual) {
		    DE fprintf (OUTPUT_FILE, "We got to 49\n");
		    // Do the output
		    outputReadLines (fwdReadName, fwdSrsData, stdout);
		    outputReadLines (revReadName, revSrsData, stdout);
		    wasOutput[fwdReadName] = wasOutput[revReadName] = 1;
		    ++identicalMateSuperReads;
		    continue;
	       }
	  }
	  DE fprintf (OUTPUT_FILE, "We got to 50\n");

	  if (numFwdFields > numRevFields) {
	       DE printf ("We got to 52\n");
	       isSubstring = 1; // To make the compiler happy
	       for (j=0; j<=numFwdFields-numRevFields; j++) {
		    DE printf ("j = %d\n", j);
		    isSubstring = 1;
		    for (k=0; 1; k++) {
			 DE printf ("k = %d, reversedRevKUnitig[k] = %d, fwdSrsDataKUnitig = %d\n", k, reversedRevSrsData.kUnitigList[k], fwdSrsData.kUnitigList[k+j]);
			 if (reversedRevSrsData.kUnitigList[k] != fwdSrsData.kUnitigList[k+j]) {
			      isSubstring = 0;
			      break; }
			 DE printf ("We got to 55\n");
			 ucChar2 = reversedRevSrsData.offsetAndOriList[k];
			 DE { outputReadLines (readName, revSrsData, stdout); }
			 DE { outputReadLines (readName, reversedRevSrsData, stdout); }
			 tOri = ucChar2 % 2;
			 DE printf ("tOri = %d\n", tOri);
			 shiftCodeAndIfAtEnd(ucChar2) {
			      tOri2 = fwdSrsData.offsetAndOriList[k+j] % 2;
			      DE printf ("tOri2 = %d\n", tOri2);
			      if (tOri != tOri2)
				   isSubstring = 0;
			      break; }
			 else {
			      if (reversedRevSrsData.offsetAndOriList[k] != fwdSrsData.offsetAndOriList[k+j]) {
				   isSubstring = 0;
				   break; } }
		    }
		    if (isSubstring)
			 break;
	       }
	       DE printf ("isSubstring = %d\n", isSubstring);
	       if (isSubstring) {
		    DE fprintf (OUTPUT_FILE, "We got to 56\n");
		    getKUnitigsFromSuperRead (&reversedRevSrsData, kUnitigList1);
		    someKUnitigIsNotRepetitive = reportIfSomeKUnitigInArrayIsNonRepetitive (kUnitigList1);
		    outputReadLines (fwdReadName, fwdSrsData, stdout);
		    if (! someKUnitigIsNotRepetitive) {
			 outputReadLines (revReadName, revSrsData, stdout); }
		    else {
			 ++numContainedSuperReadsAt1;
			 // Figure out how to reverse fwdSrsData
			 returnReversedSuperReadData (&fwdSrsData, &reversedFwdSrsData);
			 outputReadLines (revReadName, reversedFwdSrsData, stdout); }
		    wasOutput[fwdReadName] = wasOutput[revReadName] = 1;
		    continue;
	       }
	  }
	  DE printf ("We got to 60\n");
	  DE printf ("numFwdFields = %d, numRevFields = %d\n", numFwdFields, numRevFields);
	  if (numFwdFields < numRevFields) {
	       DE printf ("At 65\n");
	       isSubstring = 1; // To make the compiler happy
	       for (j=0; j<=numRevFields-numFwdFields; j++) {
		    isSubstring = 1;
		    for (k=0; 1; k++) {
			 if (fwdSrsData.kUnitigList[k] != reversedRevSrsData.kUnitigList[k+j]) {
			      isSubstring = 0;
			      break; }
			 ucChar2 = fwdSrsData.offsetAndOriList[k];
			 getOffsetAndOriFromCode(ucChar2, ucChar2, tOri);
			 if (ucChar2 < OFFSET_MEANING_END) {
			      if (fwdSrsData.offsetAndOriList[k] != reversedRevSrsData.offsetAndOriList[k+j]) {
				   isSubstring = 0;
				   break; } }
			 else {
			      tOri2 = reversedRevSrsData.offsetAndOriList[k+j] % 2;
			      if (tOri != tOri2)
				   isSubstring = 0;
			      break; }
		    }
		    if (isSubstring)
			 break;
	       }
	       DE fprintf (OUTPUT_FILE, "We got to 70\n");
	       if (isSubstring) {
		    getKUnitigsFromSuperRead (&fwdSrsData, kUnitigList1);
		    someKUnitigIsNotRepetitive = reportIfSomeKUnitigInArrayIsNonRepetitive (kUnitigList1);
		    if (! someKUnitigIsNotRepetitive) {
			 outputReadLines (fwdReadName, fwdSrsData, stdout); }
		    else {
			 ++numContainedSuperReadsAt2;
			 outputReadLines (fwdReadName, reversedRevSrsData, stdout);
		    }
		    outputReadLines (revReadName, revSrsData, stdout);
		    wasOutput[fwdReadName] = wasOutput[revReadName] = 1;
		    continue;
	       }
	  }
	  // The only time we should get here is if neither is a substring of the other
	  DE printf ("We got to 200\n");
	  // This returns badSrsDataSentinel (a constant) if it fails
	  srsResultsPtr = calculateInsertSuperRead (&fwdSrsData, &reversedRevSrsData, &srsData2);
	  if (srsResultsPtr == badSrsDataSentinel) {
	       DE printf ("We got to 210\n");
	       int canForceMerge = 0, changeLine = 0;
	       unsigned int localKUnitig = (fwdSrsData.kUnitigList)[0];
	       if ((localKUnitig == (reversedRevSrsData.kUnitigList)[0]) && (((fwdSrsData.offsetAndOriList)[0] % 2) == ((reversedRevSrsData.offsetAndOriList)[0] % 2)) && (! isRepetitive[localKUnitig]))
		    canForceMerge = 1;
	       localKUnitig = (fwdSrsData.kUnitigList)[numFwdFields-1];
	       if ((localKUnitig == (reversedRevSrsData.kUnitigList)[numRevFields-1]) && ((fwdSrsData.offsetAndOriList)[numFwdFields-1] == (reversedRevSrsData.offsetAndOriList)[numRevFields-1]) && (! isRepetitive[localKUnitig]))
		    canForceMerge = 1;
	       forwardLen = getSuperReadLengthFromSuperReadStruct (&fwdSrsData);
	       reverseLen = getSuperReadLengthFromSuperReadStruct (&reversedRevSrsData);
	       if (canForceMerge) {
		    ++forcedMerges;
		    if (forwardLen >= reverseLen)
			 changeLine = 2;
		    else
			 changeLine = 1;
//		    ++diffAndFailToMerge; // Incorrectly here before
	       }
	       // The next 2 lines are diff from perl version and now correct
	       else
		    ++diffAndFailToMerge;
	       DE fprintf (OUTPUT_FILE, "readName = %llu : ", fwdReadName);
	       if (changeLine != 1) {
		    outputReadLines (fwdReadName, fwdSrsData, stdout); }
	       else {
		    outputReadLines (fwdReadName, reversedRevSrsData, stdout); }
	       if (changeLine != 2) {
		    outputReadLines (revReadName, revSrsData, stdout); }
	       else {
		    outputReadLines (revReadName, fwdSrsData, stdout); }
	  }
	  else {
	       ++newlyMerged;
	       returnReversedSuperReadData (srsResultsPtr, &srsData4);
	       outputReadLines (fwdReadName, (*srsResultsPtr), stdout);
	       outputReadLines (revReadName, (srsData4), stdout);
	  }
	  wasOutput[fwdReadName] = wasOutput[revReadName] = 1;
     }

     for (ullTemp=0; ullTemp < numInserts; ullTemp++) {
	  unsigned long long readName;
	  if (insertCount[ullTemp] != 1) 
	       continue;
	  readName = 2 * ullTemp + insertEndForFirstReadOfInsert[ullTemp];
	  if (wasOutput[readName])
	       continue;
	  setSrsStruct (&srsGlobalData, startIndex[ullTemp], &srsData4, 0);
	  outputReadLines (readName, (srsData4), stdout);
     }
     outfile = fopen (statsFile, "w");
     fprintf (outfile, "num identical super-reads: %llu\n", identicalMateSuperReads);
     fprintf (outfile, "num contained at 1: %llu\n", numContainedSuperReadsAt1);
     fprintf (outfile, "num contained at 2: %llu\n", numContainedSuperReadsAt2);
     fprintf (outfile, "diff and merged: %llu\n", newlyMerged);
     fprintf (outfile, "forced (dovetail) mate merges: %llu\n", forcedMerges);
     fprintf (outfile, "diff but failed to merge: %llu\n", diffAndFailToMerge);
     fclose (outfile);

     return (0);
}

void setStatsFilename (char *thisProgramName)
{
     char statsFileOutDir[1000], execPrefix[1000], tempString[1000];
     unsigned int i;

     strcpy (tempString, superReadsFile);
     strcpy (statsFileOutDir, dirname (tempString));
     strcpy (tempString, thisProgramName);
     strcpy (execPrefix, basename (tempString));
     for (i=0; i<strlen(execPrefix); i++)
	  if (execPrefix[i] == '.')
	       break;
     execPrefix[i] = '\0';
     sprintf (statsFile, "%s/%s.stats.txt", statsFileOutDir, execPrefix);
}

void processArgs (int argc, char *argv[])
{
     int fileNum, i;

     minSepNonRptKUni = 1000000;
     fileNum = 0;
     for (i=1; i<argc; i++) {
	  if (strcmp (argv[i], "-min-sep-non-rpt-k-uni") == 0) {
	       ++i;
	       minSepNonRptKUni = atoi (argv[i]);
	       continue; }
	  if (fileNum == 0)
	       strcpy (superReadsFile, argv[i]);
	  else if (fileNum == 1)
	       strcpy (numKUnitigsFile, argv[i]);
	  else if (fileNum == 2)
	       strcpy (chimericReadsFile, argv[i]);
	  else if (fileNum == 3)
	       strcpy (repetitiveKUnitigsFile, argv[i]);
	  else if (fileNum == 4)
	       strcpy (kUnitigLengthsFile, argv[i]);
	  else if (fileNum == 5)
	       strcpy (numReadsFile, argv[i]);
	  ++fileNum;
     }
}

int getFldsFromLine (char *cptrHold)
{
     int numFlds=0, state = 0;
     char *cptr;

     for (cptr=cptrHold; *cptr; cptr++) {
          if (isspace (*cptr)) { state = 0;*cptr = 0; }
          else {
               if (state == 1) continue;
               flds[numFlds] = cptr;
               ++numFlds;
               state = 1;
          }
     }
     return (numFlds);
}

int getSuperReadLengthFromSuperReadStruct (struct superReadStruct *srsPtr)
{
     int localLength, i;
     unsigned char offsetAndOriChar;

     localLength = 0;
     for (i=0; 1; i++) {
	  localLength += kUnitigLengths[(srsPtr->kUnitigList)[i]];
	  offsetAndOriChar = (srsPtr->offsetAndOriList)[i];
	  shiftCodeAndIfAtEnd(offsetAndOriChar)
	       break;
	  localLength -= offsetAndOriChar;
     }

     return (localLength);
}

int getNumKUnitigsInSuperReadStruct (struct superReadStruct *srsPtr)
{
     int i;
     unsigned char offsetAndOriChar;

     for (i=0; 1; i++) {
          offsetAndOriChar = (srsPtr->offsetAndOriList)[i];
          shiftCodeAndIfAtEnd(offsetAndOriChar)
               break;
     }
     return (i+1);
}

// The following was checked inside the function
void returnReversedSuperReadData (struct superReadStruct *inputSrsPtr, struct superReadStruct *outputSrsPtr)
{
     int i, j, k;
     unsigned char ucVal, ucVal2;
     // printf ("At rR10\n");
//     printf ("Inside returnReversedSuperReadData\n");
//     printf ("Before: "); printSrsStruct (outline, inputSrsPtr); fputs (outline, stdout); fputc ('\n', stdout);
     for (i=0; 1; i++) {
	  ucVal = (inputSrsPtr->offsetAndOriList)[i];
	  // printf ("i = %d, ucVal = %d\n", i, ucVal);
	  shiftCodeAndIfAtEnd(ucVal)
	       break;
     }
     // printf ("At rR20\n");
     for (j=i, k=0; j>=0; j--, k++) {
	  // fprintf (OUTPUT_FILE, "j = %d, k = %d\n", j, k);
	  (outputSrsPtr->kUnitigList)[k] = (inputSrsPtr->kUnitigList)[j];
	  if (j>0)
	       ucVal = inputSrsPtr->offsetAndOriList[j-1] >> 1;
	  else
	       ucVal = OFFSET_MEANING_END;
	  ucVal <<= 1;
	  ucVal2 = (inputSrsPtr->offsetAndOriList)[j] % 2;
	  ucVal2 = 1 - ucVal2;
	  ucVal += ucVal2;
	  (outputSrsPtr->offsetAndOriList)[k] = ucVal;
     }
//     printf ("After: "); printSrsStruct (outline, outputSrsPtr); fputs (outline, stdout); fputc ('\n', stdout);

     return;
}

struct superReadStruct *getMinOfSuperReadAndReverseComplement (struct superReadStruct *inputSrsPtr)
{
     int i;
     unsigned int kUnitig1, kUnitig2;
     unsigned char overlap1, overlap2, ori1, ori2;

     returnReversedSuperReadData (inputSrsPtr, &srsData1);
     for (i=0; 1; i++) {
	  kUnitig1 = (inputSrsPtr->kUnitigList)[i];
	  kUnitig2 = (srsData1.kUnitigList)[i];
	  if (kUnitig1 < kUnitig2)
	       return (inputSrsPtr);
	  if (kUnitig1 > kUnitig2)
	       return (&srsData1);
	  overlap1 = (inputSrsPtr->offsetAndOriList)[i];
	  overlap2 = (srsData1.offsetAndOriList)[i];
	  getOffsetAndOriFromCode(overlap1, overlap1, ori1);
	  getOffsetAndOriFromCode(overlap2, overlap2, ori2);
	  if (ori1 < ori2)
	       return (inputSrsPtr);
	  if (ori1 > ori2)
	       return (&srsData1);
	  if (overlap1 >= OFFSET_MEANING_END)
	       break;
     }
     return (inputSrsPtr);
}

// kUnitigList in next function expects a sentinel at end >= numKUnitigs
int reportIfSomeKUnitigInArrayIsNonRepetitive (unsigned int *kUnitigList)
{
     int good, i;
     unsigned int kUnitig;
//     printf ("In reportIfSomeKUnitigInArrayIsNonRepetitive\n");

     good = 0;
     for (i=0; 1; i++) {
	  kUnitig = kUnitigList[i];
//	  printf ("kUnitig = %d\n", kUnitig);
	  if (kUnitig >= numKUnitigs)
	       break;
	  if (isRepetitive[kUnitig])
	       continue;
	  good = 1;
	  break;
     }
     return (good);
}

// The following assumes the appropriate memory has been allocated.
// The kUnitigList ends in a sentinel >= numKUnitigs
void getKUnitigsFromSuperRead (struct superReadStruct *srsPtr, unsigned int *kUnitigList)
{
     int i;
     unsigned char offsetAndOriChar;

     for (i=0; 1; i++) {
	  kUnitigList[i] = (srsPtr->kUnitigList)[i];
	  offsetAndOriChar = (srsPtr->offsetAndOriList)[i];
	  shiftCodeAndIfAtEnd(offsetAndOriChar)
	       break;
     }
     kUnitigList[i+1] = numKUnitigs+1;
}

int setSrsStruct (struct superReadStruct *srsInput, unsigned long long indexIn, struct superReadStruct *srsOutput, unsigned long long indexOut)
{
     unsigned long long indexLocal;
     unsigned char ucChar1;
     for (indexLocal=0; 1; indexLocal++) {
//	  printf ("At 500\n");
	  (srsOutput->kUnitigList)[indexOut+indexLocal] = (srsInput->kUnitigList)[indexIn+indexLocal];
	  ucChar1 = (srsInput->offsetAndOriList)[indexIn+indexLocal];
//	  printf ("ucChar1 = %d, indexLocal = %llu, indexIn = %llu, indexOut = %llu\n", ucChar1, indexLocal, indexIn, indexOut);
	  (srsOutput->offsetAndOriList)[indexOut+indexLocal] = ucChar1;
	  shiftCodeAndIfAtEnd(ucChar1)
	       break;
     }
     return (indexLocal+1);
}

void printSrsStruct (char *str, struct superReadStruct *srsPtr)
{
     char *cptr, oriChar;
     int i, numWritten;
     unsigned char overlap, ori;

     cptr = str;
     for (i=0; 1; i++) {
	  overlap = (srsPtr->offsetAndOriList)[i];
	  getOffsetAndOriFromCode(overlap, overlap, ori);
	  oriChar = ((ori == 0) ? 'F' : 'R');
	  numWritten = sprintf (cptr, "%u %c", (srsPtr->kUnitigList)[i], oriChar);
	  cptr += numWritten;
	  if (overlap >= OFFSET_MEANING_END)
	       break;
	  numWritten = sprintf (cptr, " %u ", overlap);
	  cptr += numWritten;
     }
}

struct superReadStruct *calculateInsertSuperRead (struct superReadStruct *inputSrsPtr1, struct superReadStruct *inputSrsPtr2, struct superReadStruct *outputSrsPtr)
{
     int numFlds1, numFlds2;
     unsigned int lastKUnitig1;
     int matchingIndex, i, k1, k2, k3;
     int someKUnitigIsNotRepetitive;
     unsigned char lastOri1;
     DE fprintf (OUTPUT_FILE, "In calculateInsertSuperRead\n");   //
     numFlds1 = getNumKUnitigsInSuperReadStruct (inputSrsPtr1);
     lastKUnitig1 = (inputSrsPtr1->kUnitigList)[numFlds1-1];
     lastOri1 = (inputSrsPtr1->offsetAndOriList)[numFlds1-1] % 2;
     DE fprintf (OUTPUT_FILE, "numFlds1 = %d, lastKUnitig1 = %d, lastOri1 = %d", numFlds1, lastKUnitig1, lastOri1);  //
     numFlds2 = getNumKUnitigsInSuperReadStruct (inputSrsPtr2);
     matchingIndex = -1;
     for (i=0; i<numFlds2; i++) {
	  if ((inputSrsPtr2->kUnitigList)[i] == lastKUnitig1) {
	       if ((inputSrsPtr2->offsetAndOriList)[i] % 2 == lastOri1) {
		    matchingIndex = i;
		    break; } } }
     DE fprintf (OUTPUT_FILE, "matchingIndex = %d\n", matchingIndex); //
     if (matchingIndex < 0)
	  return (badSrsDataSentinel);
     DE fprintf (OUTPUT_FILE, "Got to 100\n"); //
     for (k1=numFlds1-1, k2=matchingIndex, k3=0; (k1 >= 0) && (k2 >= 0); k1--, k2--, k3++) {
	  DE fprintf (OUTPUT_FILE, "Got to 105\n"); //
	  if ((inputSrsPtr1->kUnitigList)[k1] != (inputSrsPtr2->kUnitigList)[k2])
	       return (badSrsDataSentinel);
	  DE fprintf (OUTPUT_FILE, "Got to 107\n"); //
	  kUnitigList1[k3] = (inputSrsPtr1->kUnitigList)[k1];
	  if (((inputSrsPtr1->offsetAndOriList)[k1] != (inputSrsPtr2->offsetAndOriList)[k2]) && (k3 != 0))
	       return (badSrsDataSentinel);
     }
     DE fprintf (OUTPUT_FILE, "Got to 110, k3 = %d\n", k3); //
     kUnitigList1[k3] = numKUnitigs+1;
     DE
	  for (i=0; i<k3; i++)  //
	       fprintf (OUTPUT_FILE, "i = %d, kUnitig = %d\n", i, kUnitigList1[i]); //
     someKUnitigIsNotRepetitive = reportIfSomeKUnitigInArrayIsNonRepetitive (kUnitigList1);
     DE fprintf (OUTPUT_FILE, "someKUnitigIsNotRepetitive = %d\n", someKUnitigIsNotRepetitive); //
     if (! someKUnitigIsNotRepetitive)
	  return (badSrsDataSentinel);
     DE fprintf (OUTPUT_FILE, "Got to 120\n"); //
     k3 = 0;
     for (i=0; i<numFlds1-1; i++) {
	  (outputSrsPtr->kUnitigList)[k3] = (inputSrsPtr1->kUnitigList)[i];
	  (outputSrsPtr->offsetAndOriList)[k3] = (inputSrsPtr1->offsetAndOriList)[i];
	  ++k3;
     }
     for (i=matchingIndex; i<numFlds2; i++) {
	  (outputSrsPtr->kUnitigList)[k3] = (inputSrsPtr2->kUnitigList)[i];
	  (outputSrsPtr->offsetAndOriList)[k3] = (inputSrsPtr2->offsetAndOriList)[i];
	  ++k3;
     }
     return (outputSrsPtr);
}

int superReadIsReversedForNormalForm (struct superReadStruct *srsPtr)
{
     int numFieldsLocal, i, j, itmp1, itmpI, itmpJ;
//     char outline[100000];
//     printf ("In superReadIsReversedForNormalForm\n");
//     printf ("InputSrs: "); printSrsStruct (outline, srsPtr); fputs (outline, stdout); fputc ('\n', stdout);     
     numFieldsLocal = getNumKUnitigsInSuperReadStruct (srsPtr);
//     printf ("numFieldsLocal = %d\n", numFieldsLocal);
     for (i=0, j=numFieldsLocal-1; i<=j; i++, j--) {
	  if ((srsPtr->kUnitigList)[i] < (srsPtr->kUnitigList)[j]) {
//	       printf ("Return val1 = 0, kUniI = %d, kUniJ = %d\n", (srsPtr->kUnitigList)[i], (srsPtr->kUnitigList)[j]);
	       return (0); }
	  if ((srsPtr->kUnitigList)[i] > (srsPtr->kUnitigList)[j]) {
//	       printf ("Return val2 = 1, kUniI = %d, kUniJ = %d\n", (srsPtr->kUnitigList)[i], (srsPtr->kUnitigList)[j]);
	       return (1); }
	  itmp1 = ((srsPtr->offsetAndOriList)[i] % 2) + ((srsPtr->offsetAndOriList)[j] % 2);
	  if (itmp1 == 0) {
//	       printf ("Return val3 = 0\n");
	       return (0); }
	  if (itmp1 == 2) {
//	       printf ("Return val4 = 1\n");
	       return (1); }
	  if (i < j-1) {
	       itmpI = (srsPtr->offsetAndOriList)[i] >> 1;
	       itmpJ = (srsPtr->offsetAndOriList)[j] >> 1;
	       if (itmpI < itmpJ) {
//		    printf ("Return val11 = 0\n");
		    return (0); }
	       if (itmpI > itmpJ) {
//		    printf ("Return val12 = 1\n");
		    return (1); } }
     }
     return (0);
}

void returnSuperReadInNormalForm (struct superReadStruct *inputSrsPtr, struct superReadStruct *outputSrsPtr)
{
     int needsToBeReversed;
//     char outline[100000];

//     printf ("In returnSuperReadInNormalForm\n");
//     printf ("InputSrs: "); printSrsStruct (outline, inputSrsPtr); fputs (outline, stdout); fputc ('\n', stdout);     
     needsToBeReversed = superReadIsReversedForNormalForm (inputSrsPtr);
//     printf ("Needs to be reversed = %d\n", needsToBeReversed);
     if (! needsToBeReversed)
	  setSrsStruct (inputSrsPtr, 0, outputSrsPtr, 0);
     else
	  returnReversedSuperReadData (inputSrsPtr, outputSrsPtr);
}

