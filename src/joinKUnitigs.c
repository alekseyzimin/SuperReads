// For this exec we are using the unitig numbers starting from 0
// 2 (optional) args:
// 1) The overlaps file
// 2) The mate num info file
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "priorityHeap.h"
#include "redBlackTreesInsertOnly.h"

#define NUM_UNITIGS_FILE "numKUnitigs.txt"
#define UNITIG_LENGTHS_FILE "kUnitigLengths.txt"
#define OVERLAPS_FILE "tvgxa2.overlaps"
#define INPUT_MATE_NUM_FILE "matePairNums.forPaths.txt"

#define DEFAULT_MAX_OFFSET_CONSIDERED_SAME 5
#define MAX_OFFSET_TO_TEST 13000

#define MAX_NODES_ALLOWED 20000

#define BAD_OVL 0
#define GOOD_OVL 1
#define MAXIMAL_OVL 2
#define NEWLY_ADDED_OVL 3
#define UTG1_IN_UTG2_OVL 4
#define UTG2_IN_UTG1_OVL 5
#define UTG1_AND_UTG2_ENDS_EQUAL_OVL (UTG1_IN_UTG2_OVL + UTG2_IN_UTG1_OVL)

#define NEWLY_ADDED_MAXIMAL_UNITIG 2

#define FRONT_END 1
#define BACK_END 2

struct overlapDataStruct
{
    int unitig1;
    int unitig2;
    int ahg;
    int bhg;
    unsigned short err;		// In tenths of a percent
    char ori;
     unsigned char overlapCode;
};

struct unitigLocStruct
{
    int unitig2;
    int frontEdgeOffset;
    char ori;
};

// In the following, lastOverlappingOffset gives us the offset for the node
//  where we have the largest overlap into the node
// It is set artificially to the length of the first unitig to start
struct abbrevUnitigLocStruct
{
    int frontEdgeOffset;
    int lastOffsetAtOverlap;
    unsigned short pathNum; 
    char ori;
};

struct unitigPathPrintStruct
{
     int unitig1;
     int frontEdgeOffset;
     int numOverlapsIn;
     int numOverlapsOut;
     char ori;
};

int *startOverlapByUnitig, *unitigLengths;
unsigned char *isUnitigMaximal;
struct overlapDataStruct *overlapData;
struct unitigLocStruct *unitigLocData1, *unitigLocData2;
int maxOffsetToConsiderTheSame;
int *treeReinitList, numTreesUsed;
int *startOverlapIndexByUnitig2, *unitig2OverlapIndex;
int mateUnitig1, mateUnitig2;
unsigned char mateUnitig1ori, mateUnitig2ori;
struct heapStruct priorityQueue;
int curPathNum;
int treeSize;
int maxAllowedErrorRate;	// In tenths
int minOverlapLength;
int maxDiffInsertSizesForPrinting;
FILE *outfile;
// If the following is set then we allow overlaps through which
// fail our path-creation criteria when they cover the end of a
// unitig which would otherwise be uncovered
int expandOverlapSet, reportPaths;
char *flds[100];

FILE *Fopen (const char *fn, const char *mode);
FILE *Popen (const char *fn, const char *mode);
int getFldsFromLine (char *cptrHold);
int getInt (char *fname);
int unitigLocStructCompare (struct unitigLocStruct *ptr1,
      struct unitigLocStruct *ptr2);
int unitigLocStructCompareReversed (struct unitigLocStruct *ptr1,
      struct unitigLocStruct *ptr2);
// The following returns the overlap length if it is greater than the
// existing largest overlap on the end. It returns -1 if not.
int getOvlLenFromOvlIndicesPlus (int maxOvlIndex, int j, int maxOvlLen, int whichEnd);
int findOtherOverlapIndex (int ovlIndex1);
void printIfGood (struct abbrevUnitigLocStruct *ptr);
void completePathPrint (struct abbrevUnitigLocStruct *ptr);
void printPathNode (struct unitigPathPrintStruct *ptr);
void funcToGetTreeSize (void *ptr); // Adds 1 to treeSize (a global) each time

// RB tree data stuff
struct RBTreeStruct *treeArr, *treeArr2;
struct dataArrayStruct dataArr, dataArr2;
int abbrevLocStructCompForSearch (struct abbrevUnitigLocStruct *ptr1,
      struct abbrevUnitigLocStruct *ptr2);
int abbrevLocStructCompForSort (struct abbrevUnitigLocStruct *ptr1,
      struct abbrevUnitigLocStruct *ptr2);
int unitigPathPrintStructComp (struct unitigPathPrintStruct *ptr1,
      struct unitigPathPrintStruct *ptr2);

#ifndef mallocOrDie
#define mallocOrDie(name, num, type) fprintf (stderr, "Allocating %lu bytes for %s.\n", (unsigned long) ((num) * sizeof ( type )), #name); \
name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }
#endif

// #define DEBUG 814
// #define DEBUG 1541
// #define DEBUG 349
// #define DEBUG 727157
// $define DEBUG 41218
// #define DEBUG 25129
// #define DEBUG 5374

int main (int argc, char **argv)
{
    FILE *infile;
    char line[200], cmd[400], ori, *overlapsFn = OVERLAPS_FILE, *matePairsInFn =
	  INPUT_MATE_NUM_FILE;
    int unitig1, unitig2, overlapCount = 0, itemp, itempHold;
    int unitigNum, numUnitigs, ahg, bhg, firstUnitigNum = 0;
    int offset;
    int elementIndex;
    int i, j, *iptr;
    float err1, err2;
    struct unitigLocStruct unitigLocVal;
#if DEBUG
    struct unitigLocStruct *RLPtr;
    int unitigForDebugging = DEBUG;
#endif
    struct abbrevUnitigLocStruct abbrevUnitigLocVal, *abbRLPtr;
    void *vptr;
    int lastOffsetToTest = 6000;
    int overlapLength;
    int requireOverlapsToCoverUnitig = 0;
    int lastPositionOverlapped;
    int numFieldsInput = -1;
    int insertLengthMean=0, insertLengthStdev=0;
    /* The following is to take care of the following special case:
       We are going from mateUnitig1 to mateUnitig2, and the first unitig that
       overlaps mateUnitig2 (going from mateUnitig1) contains it. Then
       if we always moved forward we would never encounter mateUnitig2.
       isSpecialCase1 is used to keep track of this situation when it occurs
       It adds the offset of mateUnitig2 to the tree of dists for mateUnitig2,
       but doesn't add it to the priority queue. (This only occurs if
       mateUnitig2 has reverse orientation.)
     */
    int isSpecialCase1;
    int maxNodes;
    char mateUnitig1str[100], mateUnitig2str[100];
    int numNewlyAddedMaximalUnitigs;
    int numFilenamesSet = 0;
    int numFlds;

    maxAllowedErrorRate = 0; // In tenths
    minOverlapLength = 40;
    expandOverlapSet = 0;
    reportPaths = 0;
//    outfile = fopen ("/localraid/tri/out.gaps5", "w");
    outfile = stdout;

    maxDiffInsertSizesForPrinting = 5;
    for (i=1; i<argc; i++) {
	 if (argv[i][0] == '-') {
	      if (strcmp (argv[i], "-max-diff-insert-sizes-for-printing") == 0) {
		   ++i;
		   maxDiffInsertSizesForPrinting = atoi (argv[i]);
	      }
	      else if (strcmp (argv[i], "-report-paths") == 0)
		   reportPaths = 1;
	      else if (strcmp (argv[i], "-require-overlaps-to-cover-unitig") == 0)
		   requireOverlapsToCoverUnitig = 1;
	      else if (strcmp (argv[i], "-expand-overlap-set") == 0)
		   expandOverlapSet = 1;
	      else if (strcmp (argv[i], "-min-overlap-length") == 0) {
		   ++i;
		   minOverlapLength = atoi (argv[i]);
	      }
	      else if (strcmp (argv[i], "-max-allowed-error-rate") == 0) {
		   float ftemp;
		   ++i;
		   ftemp = atof (argv[i]);
		   if (ftemp < 0.12) // It's an error rate
			maxAllowedErrorRate = (int) (ftemp * 1000 + .5);
		   else  // It's reported in tenths
			maxAllowedErrorRate = (int) (ftemp + .5);
	      }
	      else { // We need to allow -h later; for now we just exit
		   return (-1);
	      }
	 }
	 else if (numFilenamesSet == 0) {
	      overlapsFn = argv[i];
	      ++numFilenamesSet; }
	 else if (numFilenamesSet == 1) {
	      matePairsInFn = argv[i];
	      ++numFilenamesSet; }
    }

    mateUnitig1ori = 'F'; mateUnitig2ori = 'R';
// Get the number of unitigs
    numUnitigs = getInt (NUM_UNITIGS_FILE);
    mallocOrDie (startOverlapByUnitig, numUnitigs + 1 + firstUnitigNum, int);
    mallocOrDie (startOverlapIndexByUnitig2, numUnitigs + 1 + firstUnitigNum, int);

    mallocOrDie (unitigLengths, numUnitigs + firstUnitigNum, int);
    // Here we read in the unitig lengths, allowing for either type of length
    // format
    infile = Fopen (UNITIG_LENGTHS_FILE, "r");
    fgets (line, 200, infile);
    numFlds = getFldsFromLine (line);
    rewind (infile);
    if (numFlds == 1) {
	 for (i = 0, iptr = unitigLengths + firstUnitigNum; i < numUnitigs; i++, iptr++)
	      fscanf (infile, "%d\n", iptr);
    }
    else {
	 while (fgets (line, 200, infile)) {
	      getFldsFromLine (line);
	      unitigLengths[atoi(flds[0])] = atoi (flds[1]);
	 }
    }
    fclose (infile);

    mallocOrDie (isUnitigMaximal, numUnitigs + firstUnitigNum, unsigned char);

// Set up space to keep the overlaps file
    sprintf (cmd, "zcat -f %s", overlapsFn);
    infile = Popen (cmd, "r");
    while (fgets (line, 200, infile))
    {
	sscanf (line, "%d %d ", &unitig1, &unitig2);
	if (unitig1 >= unitig2) 
	     continue;
	overlapCount++;
	startOverlapByUnitig[unitig1]++;
	startOverlapByUnitig[unitig2]++;
    }
    mallocOrDie (overlapData, 2 * overlapCount, struct overlapDataStruct);
    mallocOrDie (unitig2OverlapIndex, 2 * overlapCount, int);
    for (unitigNum = 1; unitigNum < numUnitigs + 1 + firstUnitigNum; unitigNum++)
	startOverlapByUnitig[unitigNum] += startOverlapByUnitig[unitigNum - 1];
    for (unitigNum = 0; unitigNum < numUnitigs + 1 + firstUnitigNum; unitigNum++)
	 startOverlapIndexByUnitig2[unitigNum] = startOverlapByUnitig[unitigNum];

// Set up the priority queue
    initializeEmptyHeap (priorityQueue, struct unitigLocStruct,
	  unitigLocStructCompare);

// Set up the RB trees
    initializeEmptyTrees (treeArr, numUnitigs + 1, dataArr,
	  struct abbrevUnitigLocStruct, abbrevLocStructCompForSort,
	  abbrevLocStructCompForSearch);
    mallocOrDie (treeReinitList, numUnitigs + 1, int);
    
// Set up the RB tree for the final paths
    initializeEmptyTreesWithDataSize (treeArr2, 1, dataArr2,
          struct unitigPathPrintStruct, 40000, unitigPathPrintStructComp,
	  unitigPathPrintStructComp);

// Unitig in the overlaps file
    pclose (infile);
    infile = Popen (cmd, "r");
    while (fgets (line, 200, infile))
    {
	sscanf (line, "%d %d %c %d %d %f %f ", &unitig1, &unitig2, &ori, &ahg, &bhg,
		&err1, &err2);
	if (unitig1 >= unitig2)
	     continue;
	startOverlapByUnitig[unitig1]--;
	itemp = startOverlapByUnitig[unitig1];
	overlapData[itemp].unitig1 = unitig1;
	overlapData[itemp].unitig2 = unitig2;
	overlapData[itemp].ori = ori;
	overlapData[itemp].ahg = ahg;
	overlapData[itemp].bhg = bhg;
	overlapData[itemp].err = 10 * err2 + 0.5;
	overlapData[itemp].overlapCode = GOOD_OVL;
	if (overlapData[itemp].err > maxAllowedErrorRate)
	     overlapData[itemp].overlapCode = BAD_OVL;
	unitig1 = overlapData[itemp].unitig1;
	if (overlapData[itemp].ahg >= 0)
	     overlapLength = unitigLengths[overlapData[itemp].unitig1] - overlapData[itemp].ahg;
	else
	     overlapLength = unitigLengths[overlapData[itemp].unitig2] + overlapData[itemp].ahg;
	if (overlapLength < minOverlapLength) overlapData[itemp].overlapCode = BAD_OVL;
	startOverlapIndexByUnitig2[unitig2]--;
	unitig2OverlapIndex[startOverlapIndexByUnitig2[unitig2]] = itemp;
	startOverlapByUnitig[unitig2]--;
	itempHold = itemp;
	itemp = startOverlapByUnitig[unitig2];
	overlapData[itemp].unitig1 = unitig2;
	overlapData[itemp].unitig2 = unitig1;
	overlapData[itemp].ori = ori;
	overlapData[itemp].err = 10 * err2 + 0.5;
	startOverlapIndexByUnitig2[unitig1]--;
	unitig2OverlapIndex[startOverlapIndexByUnitig2[unitig1]] = itemp;
	if (ori == 'N')
	{
	    overlapData[itemp].ahg = -ahg;
	    overlapData[itemp].bhg = -bhg;
	}
	else
	{
	    overlapData[itemp].ahg = bhg;
	    overlapData[itemp].bhg = ahg;
	}
	overlapData[itemp].overlapCode = overlapData[itempHold].overlapCode;
    }
    pclose (infile);

    if (expandOverlapSet) {
	 int frontIsCovered, backIsCovered, numOvlsAdded=0;
	 for (i=firstUnitigNum; i<firstUnitigNum+numUnitigs; i++) {
	      frontIsCovered = backIsCovered = 0;
	      for (j=startOverlapByUnitig[i]; j<startOverlapByUnitig[i+1]; j++) {
		   if (overlapData[j].overlapCode == BAD_OVL) continue;
		   if (overlapData[j].ahg < 0) frontIsCovered = 1;
		   if (overlapData[j].bhg > 0) backIsCovered = 1;
	      }
	      if (frontIsCovered && backIsCovered) continue;
	      for (j=startOverlapByUnitig[i]; j<startOverlapByUnitig[i+1]; j++) {
		   if ((! frontIsCovered) && (overlapData[j].ahg < 0)) {
			overlapLength = unitigLengths[overlapData[j].unitig2] + overlapData[j].ahg;
			if (overlapLength >= 40) overlapData[j].overlapCode = NEWLY_ADDED_OVL;
		   }
		   if ((! backIsCovered) && (overlapData[j].bhg > 0)) {
			overlapLength = unitigLengths[i] - overlapData[j].ahg;
			if (overlapLength >= 40) overlapData[j].overlapCode = NEWLY_ADDED_OVL;
		   }
	      }
	 }
	 // Now we make sure that both records corresponding to the newly
	 // added ovls are in
	 for (i=firstUnitigNum; i<firstUnitigNum+numUnitigs; i++) {
	      for (j=startOverlapByUnitig[i]; j<startOverlapByUnitig[i+1]; j++) {
		   if (overlapData[j].overlapCode != NEWLY_ADDED_OVL)
			continue;
		   itemp = findOtherOverlapIndex (j);
		   overlapData[itemp].overlapCode = GOOD_OVL;
		   overlapData[j].overlapCode = GOOD_OVL;
		   ++numOvlsAdded;
	      }
	 }
	 fprintf (stderr, "numOvlsAddedToCoverUnitigEnds: %d\n", numOvlsAdded);
    }

    // The following determines the maximal unitigs, i.e. the unitigs which
    // are not contained in any other. In cases where two unitigs begin and
    // end in the same point, the lower-numbered unitig is considered to
    // contain the higher-numbered unitig
    for (i=firstUnitigNum; i<firstUnitigNum+numUnitigs; i++) {
	 isUnitigMaximal[i] = 1;
	 for (j=startOverlapByUnitig[i]; j<startOverlapByUnitig[i+1]; j++) {
	      if (overlapData[j].ahg > 0) continue;
	      if (overlapData[j].bhg < 0) continue;
	      if (overlapData[j].overlapCode == BAD_OVL) continue;
	      // It is a containment overlap
	      if ((overlapData[j].ahg != 0) || (overlapData[j].bhg != 0)) {
		   isUnitigMaximal[i] = 0;
		   break;
	      }
	      // They match up end-to-end
	      if (overlapData[j].unitig2 < j) {
		   isUnitigMaximal[i] = 0;
		   break;
	      }
	 }
    }

    numNewlyAddedMaximalUnitigs = 0;
    for (i=firstUnitigNum; i<firstUnitigNum+numUnitigs; i++) {
	 if (isUnitigMaximal[i]) {
	      isUnitigMaximal[i] = NEWLY_ADDED_MAXIMAL_UNITIG;
	      ++numNewlyAddedMaximalUnitigs;
	 }
    }

    // Now we find the maximal overlaps of the maximal unitigs
    while (numNewlyAddedMaximalUnitigs > 0) {
	 numNewlyAddedMaximalUnitigs = 0;
	 for (i=firstUnitigNum; i<firstUnitigNum+numUnitigs; i++) {
	      int maxFrontOvlLen=0, maxBackOvlLen=0;
	      int maxFrontOvlIndex=0, maxBackOvlIndex=0;
	      int otherUnitig;
	      if (isUnitigMaximal[i] != NEWLY_ADDED_MAXIMAL_UNITIG) continue;
	      isUnitigMaximal[i] = 1;
	      for (j=startOverlapByUnitig[i]; j<startOverlapByUnitig[i+1]; j++) {
		   // In the following case the overlap is at the end
		   if ((overlapData[j].ahg > 0) && (overlapData[j].bhg > 0)) {
			overlapLength = getOvlLenFromOvlIndicesPlus (maxBackOvlIndex, j, maxBackOvlLen, BACK_END);
			if (overlapLength > 0) {
			     maxBackOvlLen = overlapLength;
			     maxBackOvlIndex = j;
			}
		   }
		   // In the following case the overlap is at the beginning
		   if ((overlapData[j].ahg < 0) && (overlapData[j].bhg < 0)) {
			overlapLength = getOvlLenFromOvlIndicesPlus (maxFrontOvlIndex, j, maxFrontOvlLen, FRONT_END);
			if (overlapLength > 0) {
			     maxFrontOvlLen = overlapLength;
			     maxFrontOvlIndex = j;
			}
		   }
	      }
	      if (maxFrontOvlLen > 0) {
		   overlapData[maxFrontOvlIndex].overlapCode = MAXIMAL_OVL;
		   otherUnitig = overlapData[maxFrontOvlIndex].unitig2;
		   if (! isUnitigMaximal[otherUnitig]) {
			isUnitigMaximal[otherUnitig] = NEWLY_ADDED_MAXIMAL_UNITIG;
			++numNewlyAddedMaximalUnitigs;
		   }
		   // Now we find the other record for the overlap
		   itemp = findOtherOverlapIndex (maxFrontOvlIndex);
		   overlapData[itemp].overlapCode = MAXIMAL_OVL;
	      }
	      if (maxBackOvlLen > 0) {
		   overlapData[maxBackOvlIndex].overlapCode = MAXIMAL_OVL;
		   otherUnitig = overlapData[maxBackOvlIndex].unitig2;
		   if (! isUnitigMaximal[otherUnitig]) {
			isUnitigMaximal[otherUnitig] = NEWLY_ADDED_MAXIMAL_UNITIG;
			++numNewlyAddedMaximalUnitigs;
		   }
		   // Now we find the other record for the overlap
		   itemp = findOtherOverlapIndex (maxBackOvlIndex);
		   overlapData[itemp].overlapCode = MAXIMAL_OVL;
	      }
	      // We pass here once for each unitig in each pass
	 }
    }

    // Now calculate the containment overlaps
    for (i=startOverlapByUnitig[firstUnitigNum]; i<startOverlapByUnitig[firstUnitigNum+numUnitigs]; i++) {
	 if (overlapData[i].overlapCode != GOOD_OVL) continue;
	 if ((overlapData[i].ahg<=0) && (overlapData[i].bhg>=0))
	      overlapData[i].overlapCode = UTG1_IN_UTG2_OVL;
	 if ((overlapData[i].ahg>=0) && (overlapData[i].bhg<=0)) {
	      if (overlapData[i].overlapCode == GOOD_OVL)
		   overlapData[i].overlapCode = UTG2_IN_UTG1_OVL;
	      else
		   overlapData[i].overlapCode = UTG1_AND_UTG2_ENDS_EQUAL_OVL;
	 }
    }

    infile = Fopen (matePairsInFn, "r");
    while (fgets (line, 200, infile))
    {
	if (numFieldsInput < 0)
	{
	    int status = 0;	// In a space
	    char *cptr;
	    numFieldsInput = 0;
	    for (cptr = line; *cptr; cptr++)
	    {
		if (isspace (*cptr))
		    status = 0;
		else
		{
		    if (status == 0)
			numFieldsInput++;
		    status = 1;
		}
	    }
	}
	if (numFieldsInput == 2)
	{
	    sscanf (line, "%d %d\n", &mateUnitig1, &mateUnitig2);
#ifndef NO_OUTPUT
	    fprintf (outfile, "MATE PAIR: %d %d\n", mateUnitig1, mateUnitig2);
#endif
	}
	else
	{  
	     if (numFieldsInput == 3) {
		  sscanf (line, "%d %d %d\n", &mateUnitig1, &mateUnitig2, &insertLengthMean);
		  if (insertLengthMean == 4000) {
		       insertLengthMean = 3780;
		       insertLengthStdev = 450;
		  }
		  else if (insertLengthMean == 5500) {
		       insertLengthMean = 5273;
		       insertLengthStdev = 569;
		  }
		  else {
		       insertLengthMean = 34599; // ??
		       insertLengthStdev = 3337; // ??
		  }
	     }
	     else { // numFieldsInput should equal 4
		  int len;
		  sscanf (line, "%s %s %d %d\n", mateUnitig1str, mateUnitig2str, &insertLengthMean, &insertLengthStdev);
		  mateUnitig1ori = 'F';
		  len = strlen (mateUnitig1str);
		  if (isdigit (mateUnitig1str[len-1]))
		       mateUnitig1 = atoi (mateUnitig1str);
		  else {
		       mateUnitig1ori = mateUnitig1str[len-1];
		       mateUnitig1str[len-1] = 0;
		       mateUnitig1 = atoi (mateUnitig1str);
		       mateUnitig1str[len-1] = mateUnitig1ori;
		  }
		  mateUnitig2ori = 'R';
		  len = strlen (mateUnitig2str);
		  if (isdigit (mateUnitig2str[len-1]))
		       mateUnitig2 = atoi (mateUnitig2str);
		  else {
		       mateUnitig2ori = mateUnitig2str[len-1];
		       mateUnitig2str[len-1] = 0;
		       mateUnitig2 = atoi (mateUnitig2str);
		       mateUnitig2str[len-1] = mateUnitig2ori;
		  }
		  
	     }
#ifndef NO_OUTPUT
	     fprintf (outfile, "MATE PAIR: %s %s %d %d", mateUnitig1str, mateUnitig2str, insertLengthMean, insertLengthStdev);
#endif
	     // If the unitigs are "misoriented" then we need to adjust the
	     // insertLengthMean, since we calculate as though they both had
	     // the appropriate orientation and afterwards adjust by the
	     // unitig length when the unitig start is on the other side.
	     // This is necessary since otherwise we may be more than 4
	     // standard deviations away from the mean
	     if (mateUnitig1ori != 'F')
		  insertLengthMean += unitigLengths[mateUnitig1];
	     if (mateUnitig2ori != 'R')
		  insertLengthMean += unitigLengths[mateUnitig2];
	     if (insertLengthMean > MAX_OFFSET_TO_TEST) {
		  // It's too long so we're skipping it
#ifndef NO_OUTPUT
		  fprintf (outfile, "; avg insert length too long. Skipping.\n");
#endif
		  continue;
	     }
#ifndef NO_OUTPUT
	     else 
		  fprintf (outfile, "\n");
#endif
	     lastOffsetToTest = insertLengthMean+4*insertLengthStdev;
	     if (lastOffsetToTest > MAX_OFFSET_TO_TEST)
		  lastOffsetToTest = MAX_OFFSET_TO_TEST;
	}
#if DEBUG
	if (argc > 3) unitigForDebugging = atoi (argv[3]);
	if ((mateUnitig1 != unitigForDebugging) && (mateUnitig2 != unitigForDebugging))
	     continue;
#endif
	// Adjust overlaps for mateUnitig1 if ori not 'F' and for
	//    mateUnitig2 if ori not 'R' to what they would be if they had
	//    the desired orientation. Note that mateUnitig1 is only a
	//    source and mateUnitig2 is only a sink
	if (mateUnitig1ori != 'F') { 
	     for (j=startOverlapByUnitig[mateUnitig1]; j<startOverlapByUnitig[mateUnitig1+1]; j++) {
		  int itemp=overlapData[j].ahg;
		  overlapData[j].ahg = - overlapData[j].bhg;
		  overlapData[j].bhg = - itemp;
		  if (overlapData[j].ori == 'N') overlapData[j].ori = 'I';
		  else overlapData[j].ori = 'N';
	     }
	}
	if (mateUnitig2ori != 'R') {
	     for (j=startOverlapIndexByUnitig2[mateUnitig2]; j<startOverlapIndexByUnitig2[mateUnitig2+1]; j++) {
		  int index = unitig2OverlapIndex[j];
		  if (overlapData[index].ori == 'N') overlapData[index].ori = 'I';
		  else overlapData[index].ori = 'N';
	     }
	}
	unitigLocVal.unitig2 = mateUnitig1;
	unitigLocVal.frontEdgeOffset = unitigLengths[mateUnitig1];
	unitigLocVal.ori = 'F';
	numTreesUsed = 0;
	abbrevUnitigLocVal.frontEdgeOffset = unitigLocVal.frontEdgeOffset;
	abbrevUnitigLocVal.ori = unitigLocVal.ori;
	// We may want to change the following
	abbrevUnitigLocVal.lastOffsetAtOverlap =
	      abbrevUnitigLocVal.frontEdgeOffset;
	if (treeArr[mateUnitig1].root == TREE_NIL)
	{
	    treeReinitList[numTreesUsed] = mateUnitig1;
	    numTreesUsed++;
	}
	abbrevUnitigLocVal.pathNum = 0;
	RBTreeInsertElement (treeArr + mateUnitig1, &abbrevUnitigLocVal);

	heapInsert (&priorityQueue, &unitigLocVal);
	maxNodes = 1;
	while (!heapIsEmpty (&priorityQueue))
	{
	     if (priorityQueue.heapSize > maxNodes)
		  maxNodes = priorityQueue.heapSize;
#if DEBUG
	     printf ("Starting new offset: "); fflush (stdout);
	    for (j = 0; j < priorityQueue.heapSize; j++)
	    {
		setHeapValPtr (vptr, &priorityQueue, j);
		RLPtr = vptr;
		printf ("%d ", RLPtr->frontEdgeOffset); fflush (stdout);
	    }
	    printf ("\n"); fflush (stdout);
#endif
	    heapExtractRoot (&priorityQueue, &unitigLocVal);
	    unitig1 = unitigLocVal.unitig2;
#if DEBUG
	    printf ("unitig1 = %d, ", unitig1); fflush (stdout);
#endif
	    ori = unitigLocVal.ori;
	    offset = unitigLocVal.frontEdgeOffset;
	    abbrevUnitigLocVal.frontEdgeOffset = offset;
#if DEBUG
	    printf ("offset = %d, ", offset); fflush (stdout);
#endif
	    abbrevUnitigLocVal.ori = ori;
	    elementIndex =
		  treeFindElement (treeArr + unitig1, &abbrevUnitigLocVal);
	    assert (elementIndex != TREE_NIL);
	    setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	    abbRLPtr = vptr;
	    lastPositionOverlapped = abbRLPtr->lastOffsetAtOverlap;
#if DEBUG
	    printf ("offset = %d, lastPosOvlppd = %d, ori = %c, ", offset,
		  lastPositionOverlapped, ori);  fflush (stdout);
#endif
	    for (j = startOverlapByUnitig[unitig1];
		  j < startOverlapByUnitig[unitig1 + 1]; j++)
	    {
#if DEBUG
		printf ("Starting an overlap: err = %d", overlapData[j].err);
#endif
		if (overlapData[j].overlapCode == BAD_OVL)
		     continue;
#if 0
		if (overlapData[j].err > maxAllowedErrorRate) {
#if DEBUG
		     printf (" Err rate too high\n");
#endif
		     continue;
		}
#endif
		unitig2 = overlapData[j].unitig2;
		if (overlapData[j].ahg >= 0)
		    overlapLength = unitigLengths[unitig1] - overlapData[j].ahg;
		else
		    overlapLength = unitigLengths[unitig2] + overlapData[j].ahg;
		if (overlapLength > unitigLengths[unitig1])
		     overlapLength = unitigLengths[unitig1];
		if (overlapLength > unitigLengths[unitig2])
		     overlapLength = unitigLengths[unitig2];
#if DEBUG
		printf ("; ovl len = %d", overlapLength);
#endif
#if 0		
		if (overlapLength < minOverlapLength) {
#if DEBUG
		     printf (" Ovl too short\n");
#endif
		    continue;
		}
#endif
#if DEBUG
		printf ("; unitig2 = %d, ori = %c, ahg = %d, bhg = %d\n", unitig2,
		      overlapData[j].ori, overlapData[j].ahg,
		      overlapData[j].bhg);
#endif
		isSpecialCase1 = 0;
		if (ori == 'F')
		{
		     if (overlapData[j].bhg <= 0) 
		     {
//			  if (unitig2 == mateUnitig2)
//			       isSpecialCase1 = 1;
//			  else
			       continue;
		     }
		    bhg = overlapData[j].bhg;
		    if (overlapData[j].ori == 'N')
			abbrevUnitigLocVal.ori = 'F';
		    else
			abbrevUnitigLocVal.ori = 'R';
		    abbrevUnitigLocVal.frontEdgeOffset = offset + bhg;
		}
		else
		{
		    if (overlapData[j].ahg >= 0)
		    {
//			 if (unitig2 == mateUnitig2)
//			      isSpecialCase1 = 1;
//			 else
			      continue;
		    }
		    ahg = overlapData[j].ahg;
		    if (overlapData[j].ori == 'N')
			abbrevUnitigLocVal.ori = 'R';
		    else
			abbrevUnitigLocVal.ori = 'F';
		    abbrevUnitigLocVal.frontEdgeOffset = offset - ahg;
		}
		if (isSpecialCase1 && (abbrevUnitigLocVal.ori != 'R'))
		     continue;
		if (isSpecialCase1 &&
		    (abbrevUnitigLocVal.frontEdgeOffset < lastPositionOverlapped) &&
		    (lastPositionOverlapped != unitigLengths[mateUnitig1]))
		     continue;
		abbrevUnitigLocVal.lastOffsetAtOverlap = offset;
		// Skip if the offset is too large
		if (abbrevUnitigLocVal.frontEdgeOffset >= lastOffsetToTest)
		    continue;
		if ((requireOverlapsToCoverUnitig)
		      && (abbrevUnitigLocVal.frontEdgeOffset -
			    unitigLengths[unitig2] > lastPositionOverlapped))
		    continue;
#if DEBUG
		printf ("cur front = %d\n", abbrevUnitigLocVal.frontEdgeOffset);
#endif
		// Skip if abbrevUnitigLocVal alunitigy seen for unitig2
		elementIndex =
		      treeFindElement (treeArr + unitig2, &abbrevUnitigLocVal);
		setTreeValPtr (vptr, treeArr + unitig2, elementIndex);
		abbRLPtr = vptr;
		abbRLPtr->lastOffsetAtOverlap = offset;

		if (elementIndex != TREE_NIL)
		    continue;
#if DEBUG
		printf ("Adding to the tree\n");
#endif
		// Insert this value in the priority queue
		unitigLocVal.unitig2 = unitig2;
		unitigLocVal.frontEdgeOffset = abbrevUnitigLocVal.frontEdgeOffset;
		unitigLocVal.ori = abbrevUnitigLocVal.ori;
//		if (! isSpecialCase1)
		if ((unitigLocVal.unitig2 != mateUnitig1) && (unitigLocVal.unitig2 != mateUnitig2))
		     heapInsert (&priorityQueue, &unitigLocVal);
		if (treeArr[unitig2].root == TREE_NIL)
		{
		    // If unitig2's tree never seen before
		    //    Add unitig2 to list of trees to reinit
		    treeReinitList[numTreesUsed] = unitig2;
		    numTreesUsed++;
		}
		// Add offset to list for tree
		abbrevUnitigLocVal.pathNum = 0;
		RBTreeInsertElement (treeArr + unitig2, &abbrevUnitigLocVal);
#if DEBUG
//		if ((unitig2 == mateUnitig2) && (abbrevUnitigLocVal.ori == 'R'))
		if (unitig2 == mateUnitig2)
		     printf ("Adding distance %d for the (rev oriented) mate, unitig %d\n",
			     abbrevUnitigLocVal.frontEdgeOffset, mateUnitig2);
#endif		
		//   Make sure the root of the tree is updated (if needed)
	    }			// End of going through overlaps for unitig
	    if (maxNodes > MAX_NODES_ALLOWED)
		 break;
	}			// Ends heapIsEmpty line
	// Do output for the mate pair
	//      for (j=0; j<=numUnitigs; j++)
	curPathNum = 0;
	if (treeArr[mateUnitig2].root != TREE_NIL)
	{
	     treeSize = 0;
	     inOrderTreeWalk (treeArr + mateUnitig2, treeArr[mateUnitig2].root,
			      (void (*)()) funcToGetTreeSize);
#if 0
	     printf ("treeSize = %d\n", treeSize);
#endif
#ifndef NO_OUTPUT
	    fprintf (outfile, "mateUnitig2 = %d\n", mateUnitig2);
#endif
	    fprintf (outfile, "maxNodes = %d\n", maxNodes);
	    priorityQueue.compare = (void *) unitigLocStructCompareReversed;
	    if (maxNodes <= MAX_NODES_ALLOWED)
		 inOrderTreeWalk (treeArr + mateUnitig2, treeArr[mateUnitig2].root,
//				  (void (*)(void *)) printIfGood);
				  (void (*)(void *)) completePathPrint);
	    priorityQueue.compare = (void *) unitigLocStructCompare;
	}
	// Now clean up
	// Return the overlaps to what they were if we have
	// mateUnitig1ori == 'R'  or mateUnitig2ori == 'F'
	if (mateUnitig1ori != 'F') {
	     for (j=startOverlapByUnitig[mateUnitig1]; j<startOverlapByUnitig[mateUnitig1+1]; j++) {
		  int itemp=overlapData[j].ahg;
		  overlapData[j].ahg = - overlapData[j].bhg;
		  overlapData[j].bhg = - itemp;
		  if (overlapData[j].ori == 'N') overlapData[j].ori = 'I';
		  else overlapData[j].ori = 'N';
	     }
	}
	if (mateUnitig2ori != 'R') {
	     for (j=startOverlapIndexByUnitig2[mateUnitig2]; j<startOverlapIndexByUnitig2[mateUnitig2+1]; j++) {
		  int index = unitig2OverlapIndex[j];
		  if (overlapData[index].ori == 'N') overlapData[index].ori = 'I';
		  else overlapData[index].ori = 'N';
	     }
	}
	priorityQueue.heapSize = 0;
	for (j = 0; j < numTreesUsed; j++)
	    treeArr[treeReinitList[j]].root = TREE_NIL;
	numTreesUsed = 0;
	dataArr.arraySize = 0;
	//      return (0);
    }

    if (outfile != stdout) fclose (outfile);

    return (0);
}

// The following returns the overlap length if it is greater than the
// existing largest overlap on the end. It returns -1 if not.
int getOvlLenFromOvlIndicesPlus (int maxOvlIndex, int j, int maxOvlLen, int whichEnd)
{
     int thisUnitig, otherUnitig, overlapLength;
     if (overlapData[j].overlapCode == BAD_OVL) return (-1);
     otherUnitig = overlapData[maxOvlIndex].unitig2;
     thisUnitig = overlapData[j].unitig2;
     if (whichEnd == FRONT_END)
	  overlapLength = unitigLengths[thisUnitig] - overlapData[j].ahg;
     else
	  overlapLength = unitigLengths[thisUnitig] + overlapData[j].bhg;
     if (overlapLength < maxOvlLen) return (-1);
     if (overlapLength > maxOvlLen) return (overlapLength);
     if (unitigLengths[otherUnitig] > unitigLengths[thisUnitig]) return (-1);
     if (unitigLengths[otherUnitig] < unitigLengths[thisUnitig]) return (overlapLength);
     if (otherUnitig < thisUnitig) return (-1);
     if (otherUnitig > thisUnitig) return (overlapLength);
     if (overlapData[maxOvlIndex].ori == 'N') return (-1);
     return (overlapLength);
}

int findOtherOverlapIndex (int ovlIndex1)
{
     int unitig1, unitig2, itemp;
     unitig1 = overlapData[ovlIndex1].unitig1;
     unitig2 = overlapData[ovlIndex1].unitig2;
     for (itemp=startOverlapByUnitig[unitig2]; itemp<startOverlapByUnitig[unitig2+1]; itemp++) {
	  if (overlapData[itemp].unitig2 != unitig1) continue;
	  if (overlapData[itemp].ori != overlapData[ovlIndex1].ori) continue;
	  if (overlapData[itemp].ori == 'N') {
	       if (overlapData[itemp].ahg + overlapData[ovlIndex1].ahg != 0)
		    continue;
	  }
	  else {
	       if (overlapData[itemp].ahg != overlapData[ovlIndex1].bhg)
		    continue;
	  }
	  // If we get here they agree
	  break;
     }
     return (itemp);
}

void printIfGood (struct abbrevUnitigLocStruct *ptr)
{
     int val;
     if (ptr->ori == 'R') {
	  val = ptr->frontEdgeOffset;
	  if (mateUnitig1ori != 'F') val -= unitigLengths[mateUnitig1];
	  if (mateUnitig2ori != 'R') val -= unitigLengths[mateUnitig2];
#ifndef NO_OUTPUT
	  fprintf (outfile, "%d\n", val);
#endif
     }
}

void printPathNode (struct unitigPathPrintStruct *ptr)
{
     int beginOffset, endOffset;
     if (ptr->unitig1 == mateUnitig1) ptr->ori = mateUnitig1ori;
     if (ptr->unitig1 == mateUnitig2) ptr->ori = mateUnitig2ori;
     if (ptr->ori == 'F') {
	  endOffset = ptr->frontEdgeOffset;
	  beginOffset = endOffset - unitigLengths[ptr->unitig1];
     }
     else {
	  beginOffset = ptr->frontEdgeOffset;
	  endOffset = beginOffset - unitigLengths[ptr->unitig1];
     }
     fprintf (outfile, "uni = %d, offset = %d, ori = %c, beginOffset = %d, endOffset = %d, numOvlsIn = %d, numOvlsOut = %d\n", ptr->unitig1, ptr->frontEdgeOffset, ptr->ori, beginOffset, endOffset, ptr->numOverlapsIn, ptr->numOverlapsOut);
}

void completePathPrint (struct abbrevUnitigLocStruct *ptr)
{
     struct abbrevUnitigLocStruct abbrevUnitigLocVal, *abbRLPtr;
     struct unitigLocStruct unitigLocVal;
     struct unitigPathPrintStruct unitigPathPrintVal, *rppvPtr1, *rppvPtr2;
     char ori;
     int isSpecialCase, finalOffset, minConnectingOffset, i, index;
     int unitig1, unitig2, offset, elementIndex1, elementIndex, overlapLength;
     int minConnectingUnitig=0, minConnectingOverlapIndex;
     char minConnectingOri=' ';
     void *vptr;
     // In the following we assume we move from left to right when moving from
     // mateUnitig1 to mateUnitig2.
     // First we want to find out if mateUnitig2 has an overlap which ends to
     // the left of the right-most base of mateUnitig2. If not, we put the
     // unitig which overlaps mateUnitig2 and extends the least to the left as
     // the last node of the priority queue (but we have to adjust for the 
     // missing connection between this unitig and mateUnitig2 later.)
     // Otherwise, we put mateUnitig2 on the priority queue and start.
     // Note that at this point mateUnitig2 have reverse orientation; even if
     // it was specified as forward orientation in the pair, it has been changed
     // to reverse orientation. Similarly, mateUnitig1 has forward orientation.
     ++curPathNum;
#if 0
     printf ("curPathNum = %d, nodePathNum = %d; ", curPathNum, ptr->pathNum);
#endif
     ori = ptr->ori;
     if (ori != 'R') return;
     isSpecialCase = 1;
     finalOffset = ptr->frontEdgeOffset;
     minConnectingOffset = finalOffset + 1000000;
     fprintf (outfile, "%d\n", finalOffset);
     for (i=startOverlapIndexByUnitig2[mateUnitig2]; i<startOverlapIndexByUnitig2[mateUnitig2+1]; i++) {
	  index = unitig2OverlapIndex[i];
	  if (overlapData[index].overlapCode == BAD_OVL) continue;
#if 0
	  if (overlapData[index].err > maxAllowedErrorRate)
	       continue;
#endif
	  unitig1 = overlapData[index].unitig1;
	  if (overlapData[index].ahg >= 0)
	       overlapLength = unitigLengths[unitig1] - overlapData[index].ahg;
	  else
	       overlapLength = unitigLengths[mateUnitig2] + overlapData[index].ahg;
	  if (overlapLength > unitigLengths[unitig1])
	       overlapLength = unitigLengths[unitig1];
	  if (overlapLength > unitigLengths[mateUnitig2])
	       overlapLength = unitigLengths[mateUnitig2];
#if 0
	  if (overlapLength < minOverlapLength)
	       continue;
#endif
	  if (overlapData[index].ori == 'N')
	       abbrevUnitigLocVal.ori = ori;
	  else {
	       if (ori == 'F') abbrevUnitigLocVal.ori = 'R';
	       else abbrevUnitigLocVal.ori = 'F';
	  }
	  if (abbrevUnitigLocVal.ori == 'F')
	       abbrevUnitigLocVal.frontEdgeOffset = finalOffset - overlapData[index].bhg;
	  else
	       abbrevUnitigLocVal.frontEdgeOffset = finalOffset + overlapData[index].ahg;
//	  printf ("We are at 1\n");
	  elementIndex = treeFindElement (treeArr + unitig1, &abbrevUnitigLocVal);
	  if (elementIndex == TREE_NIL) continue;
//	  printf ("We are at 2\n");
	  setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	  abbRLPtr = vptr;
#if 0
	  printf ("frontEdgeOffset = %d\n", abbRLPtr->frontEdgeOffset);
#endif
	  if (abbRLPtr->frontEdgeOffset < finalOffset) {
	       isSpecialCase = 0;
	       break;
	  }
	  else if (abbRLPtr->frontEdgeOffset < minConnectingOffset) {
	       minConnectingOffset = abbRLPtr->frontEdgeOffset;
	       minConnectingUnitig = unitig1;
	       minConnectingOri = abbrevUnitigLocVal.ori;
	       minConnectingOverlapIndex = index;
	  }
     }
     if (isSpecialCase) {
	  unitigLocVal.unitig2 = minConnectingUnitig;
	  unitigLocVal.frontEdgeOffset = minConnectingOffset;
	  unitigLocVal.ori = minConnectingOri;
	  // ..and for the path
	  unitigPathPrintVal.unitig1 = mateUnitig2;
	  unitigPathPrintVal.numOverlapsIn = 1;
	  unitigPathPrintVal.numOverlapsOut = 0;
	  unitigPathPrintVal.ori = 'R'; // Forced; may be adjusted later
	  unitigPathPrintVal.frontEdgeOffset = finalOffset;
	  RBTreeInsertElement (treeArr2, &unitigPathPrintVal);
	  unitigPathPrintVal.unitig1 = minConnectingUnitig;
	  unitigPathPrintVal.frontEdgeOffset = minConnectingOffset;
	  unitigPathPrintVal.numOverlapsIn = 0;
	  unitigPathPrintVal.numOverlapsOut = 1;
	  unitigPathPrintVal.ori = minConnectingOri;
	  RBTreeInsertElement (treeArr2, &unitigPathPrintVal);
     }
     else {
	  unitigLocVal.unitig2 = mateUnitig2;
	  unitigLocVal.frontEdgeOffset = finalOffset;
	  unitigLocVal.ori = 'R';
	  unitigPathPrintVal.unitig1 = mateUnitig2;
	  unitigPathPrintVal.frontEdgeOffset = finalOffset;
	  unitigPathPrintVal.numOverlapsIn = 0;
	  unitigPathPrintVal.numOverlapsOut = 0;
	  unitigPathPrintVal.ori = 'R';
#if 0
	  printf ("Inserting unitig1 = %d, offset = %d, ori = %c\n", unitigPathPrintVal.unitig1, unitigPathPrintVal.frontEdgeOffset, unitigPathPrintVal.ori);
#endif
	  RBTreeInsertElement (treeArr2, &unitigPathPrintVal);
     }
#if 0
     printf ("isSpecialCase = %d, unitigLocVal = %d, %d, %c\n", isSpecialCase, mateUnitig2, finalOffset, unitigLocVal.ori);
#endif
     heapInsert (&priorityQueue, &unitigLocVal);
     while (!heapIsEmpty (&priorityQueue)) {
	  heapExtractRoot (&priorityQueue, &unitigLocVal);
	  unitig2 = unitigLocVal.unitig2;
	  offset = unitigLocVal.frontEdgeOffset;
	  ori = unitigLocVal.ori;
	  unitigPathPrintVal.unitig1 = unitig2;
	  unitigPathPrintVal.frontEdgeOffset = offset;
	  unitigPathPrintVal.ori = ori;
	  elementIndex1 = treeFindElement (treeArr2, &unitigPathPrintVal);
	  setTreeValPtr (vptr, treeArr2, elementIndex1);
	  rppvPtr1 = vptr;
#if 0
	  printf ("unitig2 = %d, offset = %d, ori = %c; elementIndex1 = %d\n", unitig2, offset, ori, elementIndex1);
#endif
	  for (i=startOverlapIndexByUnitig2[unitig2]; i<startOverlapIndexByUnitig2[unitig2+1]; i++) {
	       index = unitig2OverlapIndex[i];
	       unitig1 = overlapData[index].unitig1;
	       if (overlapData[index].overlapCode == BAD_OVL)
		    continue;
#if 0
	       if (overlapData[index].err > maxAllowedErrorRate)
		    continue;
#endif
	       if (overlapData[index].ahg >= 0)
		    overlapLength = unitigLengths[unitig1] - overlapData[index].ahg;
	       else 
		    overlapLength = unitigLengths[unitig2] + overlapData[index].ahg;
	       if (overlapLength > unitigLengths[unitig1])
		    overlapLength = unitigLengths[unitig1];
	       if (overlapLength > unitigLengths[unitig2])
		    overlapLength = unitigLengths[unitig2];
#if 0
	       if (overlapLength < minOverlapLength)
		    continue;
#endif
	       if (unitig1 == mateUnitig2) continue;
	       if (overlapData[index].ori == 'N')
		    abbrevUnitigLocVal.ori = ori;
	       else {
		    if (ori == 'F')
			 abbrevUnitigLocVal.ori = 'R';
		    else
			 abbrevUnitigLocVal.ori = 'F';
	       }
	       if (abbrevUnitigLocVal.ori == 'F')
		    abbrevUnitigLocVal.frontEdgeOffset = offset - overlapData[index].bhg;
	       else
		    abbrevUnitigLocVal.frontEdgeOffset = offset + overlapData[index].ahg;
	       elementIndex = treeFindElement (treeArr + unitig1, &abbrevUnitigLocVal);
//	       printf ("Got to 21\n");
	       if (elementIndex == TREE_NIL) continue;
	       setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	       abbRLPtr = vptr;
//	       printf ("Got to 22\n");
	       if (abbRLPtr->frontEdgeOffset >= offset) continue;
	       if ((unitig1 == mateUnitig1) && (abbRLPtr->frontEdgeOffset > unitigLengths[mateUnitig1]))
		    continue;
//	       printf ("Got to 23, abbRLVpathNum = %d, curPathNum = %d\n", abbRLPtr->pathNum, curPathNum);
	       // It hasn't been seen in the retrace, so put on the queue
	       if (abbRLPtr->pathNum < curPathNum) {
#if 0
		    printf ("Adding node: unitig2 = %d, offset = %d, ori = %c\n", unitigLocVal.unitig2, unitigLocVal.frontEdgeOffset, unitigLocVal.ori);
#endif
		    abbRLPtr->pathNum = curPathNum;
		    unitigLocVal.unitig2 = unitig1;
		    unitigLocVal.frontEdgeOffset = abbRLPtr->frontEdgeOffset;
		    unitigLocVal.ori = abbRLPtr->ori;
		    heapInsert (&priorityQueue, &unitigLocVal);
		    unitigPathPrintVal.unitig1 = unitig1;
		    unitigPathPrintVal.frontEdgeOffset = unitigLocVal.frontEdgeOffset;
		    unitigPathPrintVal.ori = unitigLocVal.ori;
		    unitigPathPrintVal.numOverlapsIn = 0;
		    unitigPathPrintVal.numOverlapsOut = 0;
		    RBTreeInsertElement (treeArr2, &unitigPathPrintVal);
	       }
		    
	       unitigPathPrintVal.unitig1 = unitig1;
	       unitigPathPrintVal.frontEdgeOffset = abbRLPtr->frontEdgeOffset;
	       unitigPathPrintVal.ori = abbRLPtr->ori;
	       elementIndex = treeFindElement (treeArr2, &unitigPathPrintVal);
	       setTreeValPtr (vptr, treeArr2, elementIndex);
	       rppvPtr2 = vptr;
	       ++(rppvPtr2->numOverlapsOut);
	       // The following must be recalced in case the array had to be
	       // moved due to needing more space
	       setTreeValPtr (vptr, treeArr2, elementIndex1);
	       rppvPtr1 = vptr;
	       ++(rppvPtr1->numOverlapsIn);
//	       printf ("Got to 24\n");
	       
	  }
     }
#if 0
     printf ("tree root = %d\n", treeArr2[0].root);
#endif
     if ((treeSize <= maxDiffInsertSizesForPrinting) && (reportPaths))
	  inOrderTreeWalk (treeArr2, treeArr2[0].root,
			   (void (*)(void *)) printPathNode);
     
#if 0
     printf ("final offset = %d, arraySize = %d\n", finalOffset, dataArr2.arraySize);
#endif
     treeArr2[0].root = TREE_NIL;
     dataArr2.arraySize = 0;
}

void funcToGetTreeSize (void *ptr)
{
     ++treeSize;
}

int unitigLocStructCompare (struct unitigLocStruct *ptr1,
      struct unitigLocStruct *ptr2)
{
    if (ptr1->frontEdgeOffset > ptr2->frontEdgeOffset)
	return (-1);
    if (ptr1->frontEdgeOffset < ptr2->frontEdgeOffset)
	return (1);
    if (ptr1->unitig2 > ptr2->unitig2)
	return (-1);
    if (ptr1->unitig2 < ptr2->unitig2)
	return (1);
    if (ptr1->ori > ptr2->ori)
	return (-1);
    if (ptr1->ori < ptr2->ori)
	return (1);
    return (0);
}

int unitigLocStructCompareReversed (struct unitigLocStruct *ptr1,
      struct unitigLocStruct *ptr2)
{
    if (ptr1->frontEdgeOffset < ptr2->frontEdgeOffset)
	return (-1);
    if (ptr1->frontEdgeOffset > ptr2->frontEdgeOffset)
	return (1);
    if (ptr1->unitig2 < ptr2->unitig2)
	return (-1);
    if (ptr1->unitig2 > ptr2->unitig2)
	return (1);
    if (ptr1->ori < ptr2->ori)
	return (-1);
    if (ptr1->ori > ptr2->ori)
	return (1);
    return (0);
}

int abbrevLocStructCompForSort (struct abbrevUnitigLocStruct *ptr1,
      struct abbrevUnitigLocStruct *ptr2)
{
    if (ptr1->ori < ptr2->ori)
	return (-1);
    if (ptr1->ori > ptr2->ori)
	return (1);
    if (ptr1->frontEdgeOffset < ptr2->frontEdgeOffset)
	return (-1);
    if (ptr1->frontEdgeOffset > ptr2->frontEdgeOffset)
	return (1);
    return (0);
}

int abbrevLocStructCompForSearch (struct abbrevUnitigLocStruct *ptr1,
      struct abbrevUnitigLocStruct *ptr2)
{
    if (ptr1->ori < ptr2->ori)
	return (-1);
    if (ptr1->ori > ptr2->ori)
	return (1);
    if (ptr1->frontEdgeOffset <
	  ptr2->frontEdgeOffset - DEFAULT_MAX_OFFSET_CONSIDERED_SAME)
	return (-1);
    if (ptr1->frontEdgeOffset >
	  ptr2->frontEdgeOffset + DEFAULT_MAX_OFFSET_CONSIDERED_SAME)
	return (1);
    return (0);
}

int unitigPathPrintStructComp (struct unitigPathPrintStruct *ptr1,
      struct unitigPathPrintStruct *ptr2)
{
     if (ptr1->unitig1 == mateUnitig2) {
	  if (ptr1->unitig1 == ptr2->unitig1)
	       return (0);
	  else
	       return(1);
     }
     if (ptr2->unitig1 == mateUnitig2) return(-1);
     if (ptr1->frontEdgeOffset < ptr2->frontEdgeOffset) return (-1);
     if (ptr1->frontEdgeOffset > ptr2->frontEdgeOffset) return (1);
     if (ptr1->unitig1 < ptr2->unitig1) return (-1);
     if (ptr1->unitig1 > ptr2->unitig1) return (1);
     if (ptr1->ori < ptr2->ori) return (-1);
     if (ptr1->ori > ptr2->ori) return (1);
     return (0);
}

FILE *Fopen (const char *fn, const char *mode)
{
    FILE *result;
    result = fopen (fn, mode);
    if (result == NULL)
    {
	fprintf (stderr, "Couldn't open file '%s' for ", fn);
	switch (mode[0])
	{
	case 'r':
	    fprintf (stderr, "reading");
	    break;
	case 'w':
	    fprintf (stderr, "writing");
	    break;
	case 'a':
	    fprintf (stderr, "appending");
	    break;
	default:
	    fprintf (stderr, "unknown operation code '%c'", mode[0]);
	    break;
	}
	fprintf (stderr, ". Bye!\n");
	exit (-1);
    }

    return (result);
}

FILE *Popen (const char *fn, const char *mode)
{
    FILE *result;
    result = popen (fn, mode);
    if (result == NULL)
    {
	fprintf (stderr, "Couldn't open file '%s' for ", fn);
	switch (mode[0])
	{
	case 'r':
	    fprintf (stderr, "reading");
	    break;
	case 'w':
	    fprintf (stderr, "writing");
	    break;
	case 'a':
	    fprintf (stderr, "appending");
	    break;
	default:
	    fprintf (stderr, "unknown operation code '%c'", mode[0]);
	    break;
	}
	fprintf (stderr, ". Bye!\n");
	exit (-1);
    }

    return (result);
}

int getInt (char *fname)
{
    FILE *infile;
    int tval;

    infile = Fopen (fname, "r");
    fscanf (infile, "%d\n", &tval);
    fclose (infile);
    return (tval);
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

