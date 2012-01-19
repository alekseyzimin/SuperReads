// For usage, see --help

#define NEW_STUFF // Put in to get node-to-node connections
// #define KILLED111115
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <charb.hpp>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <set>
#include <map>
#include <utility>
#include <vector>
#include <list>
#include <stack>

#include <err.hpp>
#include <heap.hpp>
#include <exp_buffer.hpp>
#include <src2/joinKUnitigs_v3.hpp>
extern "C" {
#include <src2/redBlackTreesInsertOnly.h>
}

#define DEFAULT_MAX_OFFSET_CONSIDERED_SAME 5
#define MAX_OFFSET_TO_TEST 10000

#define MAX_NODES_ALLOWED 4000

#define NEWLY_ADDED_MAXIMAL_UNITIG 2

#define FRONT_END 1
#define BACK_END 2

struct overlapDataStruct
{
     int unitig1;
     int unitig2;
     int ahg;
     int bhg;
     char ori;
};

struct unitigLocStruct
{
     int unitig2;
     int frontEdgeOffset;
     char ori;
  bool operator<(const unitigLocStruct &rhs) const {
    if(frontEdgeOffset > rhs.frontEdgeOffset)
      return true;
    if(frontEdgeOffset < rhs.frontEdgeOffset)
      return false;
    if (unitig2 > rhs.unitig2)
      return true;
    if (unitig2 < rhs.unitig2)
      return false;
    if (ori > rhs.ori)
      return true;
    if (ori < rhs.ori)
      return false;
    return false;
  }
  bool operator>(const unitigLocStruct &rhs) const {
    return !operator<(rhs);
  }
};

// In the following, lastOverlappingOffset gives us the offset for the node
//  where we have the largest overlap into the node
// It is set artificially to the length of the first unitig to start
struct abbrevUnitigLocStruct
{
     int frontEdgeOffset;
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

struct kuniToReadMatchStruct
{
     int matchRgnBegin;
     int matchRgnEnd;
     int ahg;
     int bhg;
     int kUnitigMatchBegin;
     int kUnitigMatchEnd;
     int orientedReadMatchBegin;
     int orientedReadMatchEnd;
     int readMatchBegin;
     int readMatchEnd;
     int kUnitigLength;
     int readLength;
     int kUnitigNumber;
     char ori;
};
ExpBuffer<struct kuniToReadMatchStruct> evenReadMatchStructs, oddReadMatchStructs;
ExpandingBuffer<unsigned char> matchStructIsUsed;
struct unitigConnectionsForPathStruct
{
     int unitig1;
     int unitig2;
     int frontEdgeOffset1;
     int frontEdgeOffset2;
     char ori1;
     char ori2;
};
ExpandingBuffer<struct unitigConnectionsForPathStruct> unitigConnectionsForPathData;

struct augmentedUnitigPathPrintStruct
{
     int unitig1;
     int frontEdgeOffset;
     int numOverlapsIn;
     int numOverlapsOut;
     int beginOffset;
     int endOffset;
     char ori;
};
ExpandingBuffer<struct augmentedUnitigPathPrintStruct> augmentedUnitigPathPrintData;
int numUnitigPathPrintRecsOnPath;
ExpandingBuffer<int> fwdStartIndices, revStartIndices, fwdNumIndices, revNumIndices;
ExpandingBuffer<int> newNodeNumsFromOld;

int *startOverlapByUnitig, *unitigLengths;
unsigned char *isUnitigMaximal;
struct overlapDataStruct *overlapData;
struct unitigLocStruct *unitigLocData1, *unitigLocData2;
int maxOffsetToConsiderTheSame;
int *treeReinitList, numTreesUsed;
int *startOverlapIndexByUnitig2, *unitig2OverlapIndex;
int mateUnitig1, mateUnitig2;
unsigned char mateUnitig1ori, mateUnitig2ori;
int beginUnitig, endUnitig;
unsigned char beginUnitigOri, endUnitigOri;
typedef heap<unitigLocStruct>::min min_heap;
typedef heap<unitigLocStruct>::max max_heap;
min_heap forward_path_unitigs;
max_heap backward_path_unitigs;
std::set<unitigLocStruct> startingNodes, endingNodes;
int startingNodeNumber;
std::vector<unitigLocStruct> nodeArray;
std::map<unitigLocStruct, int> nodeToIndexMap;
std::set<std::pair<int, int> > edgeList, sortedEdgeList;
std::vector<int> unitigNodeNumbersForPath;
// std::vector<int, std::set <int> > fwdEdgeList, revEdgeList;
struct nodePair {
     int node1;
     int node2;
};
std::vector<struct nodePair> fwdConnections, revConnections;

int curPathNum;
int treeSize;
int minOverlapLength;
int maxDiffInsertSizesForPrinting;
int maxTotAllowableMissingOnEnds;
FILE *outfile, *outputFile;
char *flds[1000];
charb outputString(200), stderrOutputString(200);
double mean[256][256], stdev[256][256];
char rdPrefix[3], rdPrefixHold[3];
long long readNum, readNumHold;
int approxNumPaths;
double insertLengthMeanBetweenKUnisForInsertGlobal, insertLengthStdevGlobal;
// The following keeps track of the distance the 2 read mates are from the
// ends of the k-unitigs at the end
int lengthAdjustment1, lengthAdjustment2;
char superReadName[100000];
int splitJoinWindowMin, splitJoinWindowMax;
int numUnitigConnectionsForPathData;
int numPairsInOneUnitig, numSimplyJoinable, numJoinableAfterRead1Analysis,
     numJoinableAfterBothReadAnalysis, numJoinableUnresolvedAtEnd;

bool firstNodeSort (struct nodePair val1, struct nodePair val2);
bool secondNodeSort (struct nodePair val1, struct nodePair val2);

bool unitigLocStructCompare (struct unitigLocStruct ptr1,
			    struct unitigLocStruct ptr2);
bool unitigLocStructCompareReversed (struct unitigLocStruct ptr1,
				    struct unitigLocStruct ptr2);
void generateSuperReadPlacementLinesForJoinedMates (void);
int joinKUnitigsFromMates (int insertLengthMean, int insertLengthStdev);
FILE *Fopen (const char *fn, const char *mode);
FILE *Popen (const char *fn, const char *mode);
// The following returns the overlap length if it is greater than the
// existing largest overlap on the end. It returns -1 if not.
int getOvlLenFromOvlIndicesPlus (int maxOvlIndex, int j, int maxOvlLen, int whichEnd);
int findOtherOverlapIndex (int ovlIndex1);
void printIfGood (struct abbrevUnitigLocStruct *ptr);
void completePathPrint (struct abbrevUnitigLocStruct *ptr);
void printPathNode (struct unitigPathPrintStruct *ptr);
int setSuperReadNameFromAugmentedPath (void);
int getSuperReadLength(void);
void funcToGetTreeSize (abbrevUnitigLocStruct *ptr); // Adds 1 to treeSize (a global) each time
void findSingleReadSuperReads(char *readName);
void getSuperReadsForInsert (void);
int processKUnitigVsReadMatches (char *inputFilename, char *outputFilename);
int getInt (const char *fname);

// extern "C" {
     int getFldsFromLine (char *cptrHold);

// RB tree data stuff
     struct RBTreeStruct *treeArr, *treeArr2;
     struct dataArrayStruct dataArr, dataArr2;
     int abbrevLocStructCompForSearch (struct abbrevUnitigLocStruct *ptr1,
				       struct abbrevUnitigLocStruct *ptr2);
     int abbrevLocStructCompForSort (struct abbrevUnitigLocStruct *ptr1,
				     struct abbrevUnitigLocStruct *ptr2);
     int unitigPathPrintStructComp (struct unitigPathPrintStruct *ptr1,
				    struct unitigPathPrintStruct *ptr2);
// }

#ifndef mallocOrDie
#define mallocOrDie(name, num, type) fprintf (stderr, "Allocating %lu bytes for %s.\n", (unsigned long) ((num) * sizeof ( type )), #name); \
     name = (type *) calloc (num, sizeof ( type ));			\
     if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }
#endif

// #define DEBUG

int main (int argc, char **argv)
{
     joinKUnitigs_v3 args(argc, argv);
     FILE *infile;
     charb line(2000),readVsKUnitigFileName(256), outputFileName(256);
     int unitig1, unitig2, overlapCount = 0;
     int unitigNum, numUnitigs, firstUnitigNum = 0;
     int i, *iptr;
     int numFlds;

     maxTotAllowableMissingOnEnds = 2;
     minOverlapLength = 40;
#if KILLED111115
     outfile = stdout;
#endif

     maxDiffInsertSizesForPrinting = 5;
     minOverlapLength              = args.min_overlap_length_arg;

     // for (i=1; i<argc; i++) {
     //      if (argv[i][0] == '-') {
     //           if (strcmp (argv[i], "-max-diff-insert-sizes-for-printing") == 0) {
     //    	    ++i;
     //    	    maxDiffInsertSizesForPrinting = atoi (argv[i]);
     //           }
     //           else if (strcmp (argv[i], "-report-paths") == 0)
     //    	    reportPaths = 1;
     //           else if (strcmp (argv[i], "-min-overlap-length") == 0) {
     //    	    ++i;
     //    	    minOverlapLength = atoi (argv[i]);
     //           }
     //           else if (strcmp (argv[i], "-mean-and-stdev-by-prefix-file") == 0) {
     //    	    ++i;
     //    	    meanAndStdevByPrefixFn = argv[i];
     //           }
     //           else if (strcmp (argv[i], "-unitig-lengths-file") == 0) {
     //    	    ++i;
     //    	    unitigLengthsFn = argv[i];
     //           }
     //           else if (strcmp (argv[i], "-num-kunitigs-file") == 0) {
     //    	    ++i;
     //    	    numKUnitigsFn = argv[i];
     //           }
     //           else if (strcmp (argv[i], "-overlaps-file") == 0) {
     //    	    ++i;
     //    	    overlapsFn = argv[i];
     //           }
     //           else if (strcmp (argv[i], "-num-file-names") == 0) {
     //    	    ++i;
     //    	    numFilenames = atoi (argv[i]);
     //           }	       
     //           else if (strcmp (argv[i], "-prefix") == 0) {
     //    	    ++i;
     //    	    outputPrefix = argv[i];
     //           }
     //           else { // We need to allow -h later; for now we just exit
     //    	    fprintf (stderr, "Unrecognized flag %s. Bye.\n", argv[i]);
     //    	    return (-1);
     //           }
     //      }
     //      else{
     //           readVsKUnitigFile = argv[i];
     //           if(numFilenames==0)
     //    		numFilenames=1; }
     // }

     rdPrefix[2] = rdPrefixHold[2] = 0;
     infile = Fopen (args.mean_and_stdev_by_prefix_file_arg, "r");
     while (fgets (line, 2000, infile)) {
	  getFldsFromLine(line);
	  mean[(int)flds[0][0]][(int)flds[0][1]] = atof (flds[1]);
	  stdev[(int)flds[0][0]][(int)flds[0][1]] = atof (flds[2]);
     }
     fclose (infile);

     numPairsInOneUnitig = numSimplyJoinable = numJoinableAfterRead1Analysis =
	  numJoinableAfterBothReadAnalysis = numJoinableUnresolvedAtEnd = 0;

     mateUnitig1ori = 'F'; mateUnitig2ori = 'R';
// Get the number of unitigs
     numUnitigs = getInt (args.num_kunitigs_file_arg) + 1;
     mallocOrDie (startOverlapByUnitig, numUnitigs + 1 + firstUnitigNum, int);
     mallocOrDie (startOverlapIndexByUnitig2, numUnitigs + 1 + firstUnitigNum, int);

     mallocOrDie (unitigLengths, numUnitigs + 1 + firstUnitigNum, int);
     // Here we read in the unitig lengths, allowing for either type of length
     // format
     infile = Fopen (args.unitig_lengths_file_arg, "r");
     if (! fgets (line, 2000, infile))
       die << "File '" << args.unitig_lengths_file_arg << "' is of length 0 or can't be read";
     numFlds = getFldsFromLine (line);
     rewind (infile);
     if (numFlds == 1) {
	  int retCode; // For the stupid new compiler
	  for (i = 0, iptr = unitigLengths + firstUnitigNum; i < numUnitigs; i++, iptr++)
	       retCode = fscanf (infile, "%d\n", iptr);
     }
     else {
	  while (fgets (line, 2000, infile)) {
	       getFldsFromLine (line);
	       unitigLengths[atoi(flds[0])] = atoi (flds[1]);
	  }
     }
     fclose (infile);

     mallocOrDie (isUnitigMaximal, numUnitigs + firstUnitigNum, unsigned char);



// Set up space to keep the overlaps file
     int fd = open(args.overlaps_file_arg, O_RDONLY);
     if(fd == -1) {
	  perror("open failed");
	  exit(1);
     }
     struct stat stat_buf;
     if(fstat(fd, &stat_buf) == -1) {
	  perror("stat failed");
	  exit(1);
     }
     overlapData = (struct overlapDataStruct *)mmap(0, stat_buf.st_size, PROT_READ, MAP_SHARED, fd, 0);
     if(overlapData == MAP_FAILED) {
	  perror("mmap failed");
	  exit(1);
     }
     close(fd);

     // Force it in memory by touching every page
     char *end = (char*)overlapData + stat_buf.st_size;
     char whatever = 0;
     for(char *ptr = (char*)overlapData; ptr < end; ptr += getpagesize())
	  whatever ^= *ptr;

     overlapCount=int((double)stat_buf.st_size/(double)sizeof(struct overlapDataStruct)+.01);
   
     for(int j=0;j<overlapCount;j++)
     {
	  int unitig1=overlapData[j].unitig1;
	  int unitig2=overlapData[j].unitig2;
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
     mallocOrDie (unitig2OverlapIndex, overlapCount, int);
     for (unitigNum = 1; unitigNum < numUnitigs + 1 + firstUnitigNum; unitigNum++)
	  startOverlapByUnitig[unitigNum] += startOverlapByUnitig[unitigNum - 1];
     for (unitigNum = 0; unitigNum < numUnitigs + 1 + firstUnitigNum; unitigNum++)
	  startOverlapIndexByUnitig2[unitigNum] = startOverlapByUnitig[unitigNum];

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
     for(int j=0;j<overlapCount;j++)
     {
	  unitig1=overlapData[j].unitig1;
	  unitig2=overlapData[j].unitig2;

	  startOverlapByUnitig[unitig1]--;
	  startOverlapIndexByUnitig2[unitig2]--;
	  unitig2OverlapIndex[startOverlapIndexByUnitig2[unitig2]] = j;
     }

#if 0
     FILE *outfile = Fopen ("overlapDataArrayFile.txt", "w");
     for (int j=0; j<overlapCount; j++)
	  fprintf (outfile, "%d %d %c %d %d\n", overlapData[j].unitig1, overlapData[j].unitig2, overlapData[j].ori, overlapData[j].ahg, overlapData[j].bhg);
     fclose (outfile);
     outfile = Fopen ("unitig2OverlapIndex.txt", "w");
     for (int j=0; j<overlapCount; j++)
	  fprintf (outfile, "%d %d\n", j, unitig2OverlapIndex[j]);
     fclose (outfile);
     outfile = Fopen ("startOverlapByUnitig.txt", "w");
     for (int j=0; j<numUnitigs+1+firstUnitigNum; j++)
	  fprintf (outfile, "%d %d\n", j, startOverlapByUnitig[j]);
     fclose (outfile);
     outfile = Fopen ("startOverlapIndexByUnitig2.txt", "w");
     for (int j=0; j<numUnitigs+1+firstUnitigNum; j++)
	  fprintf (outfile, "%d %d\n", j, startOverlapIndexByUnitig2[j]);
     fclose (outfile);
#endif
     int ret;
     for(int i = 0; i < args.num_file_names_arg; ++i) {
	  switch(fork()) {
	  case -1:
	       perror("fork failed");
	       exit(1);
	       
	  case 0:
               sprintf(readVsKUnitigFileName,"%s_%d",args.input_prefix_arg,i);
               sprintf(outputFileName,"%s_%d",args.prefix_arg,i);
	       ret=processKUnitigVsReadMatches(readVsKUnitigFileName,outputFileName);
	       fprintf (stderr, "Num pairs with both reads in same unitig: %d\nNum pairs uniquely joinable: %d\nNum pairs after disambiguation to beginning of insert: %d\nNum pairs after disambiguation to end of insert: %d\nNum still joinable but not uniquely joinable: %d\n", numPairsInOneUnitig, numSimplyJoinable, numJoinableAfterRead1Analysis, numJoinableAfterBothReadAnalysis, numJoinableUnresolvedAtEnd);
	       exit(ret);
	       
	  default:
	       break;
	  }
     }
     
     //processKUnitigVsReadMatches (readVsKUnitigFile);
     
     int status;
     for(int i = 0; i < args.num_file_names_arg; ++i) {
	  if(wait(&status) == -1) {
	       perror("wait failed");
	       exit(1);
	  }
	  if(WIFEXITED(status)) {
	       fprintf(stderr,"sub %d exit status %d\n", i, WEXITSTATUS(status));
	  } else if(WIFSIGNALED(status)) {
	       fprintf(stderr,"sub %d signaled %d coredumped %d\n",
		      i, WTERMSIG(status), WCOREDUMP(status));
	  } else {
	       fprintf(stderr,"sub %d at a loss\n", i);
	  }
     }
     
     munmap(overlapData, stat_buf.st_size);
     
     return (0);
}

void updateMatchRecords(int readNum, char *cptr, char *flds[]) {
  ExpBuffer<struct kuniToReadMatchStruct> *structs;
  if (readNum % 2 == 0) 
    structs = &evenReadMatchStructs;
  else
    structs = &oddReadMatchStructs;
  structs->push_back(kuniToReadMatchStruct());
  struct kuniToReadMatchStruct &kUTRMS = structs->back();
  cptr = flds[1]+1;
  kUTRMS.matchRgnBegin = atoi (cptr);		
  kUTRMS.matchRgnEnd = atoi (flds[2]);				
  cptr = flds[3]+1; kUTRMS.ahg = atoi (cptr);			
  while(*cptr != ',') ++cptr;	++cptr;                                 
  kUTRMS.bhg = atoi (cptr);                                         
  kUTRMS.kUnitigMatchBegin = atoi (flds[4])-1;			
  kUTRMS.kUnitigMatchEnd = atoi (flds[5]);                          
  kUTRMS.orientedReadMatchBegin = atoi (flds[6]);                   
  kUTRMS.orientedReadMatchEnd = atoi (flds[7]);			
  if (kUTRMS.orientedReadMatchBegin < kUTRMS.orientedReadMatchEnd) { 
    --kUTRMS.orientedReadMatchBegin;				
    kUTRMS.ori = 'F';						
    kUTRMS.readMatchBegin = kUTRMS.orientedReadMatchBegin;      
    kUTRMS.readMatchEnd = kUTRMS.orientedReadMatchEnd; }	
  else {								
    --kUTRMS.orientedReadMatchEnd;                                  
    kUTRMS.ori = 'R';						
    kUTRMS.readMatchBegin = kUTRMS.orientedReadMatchEnd;	
    kUTRMS.readMatchEnd = kUTRMS.orientedReadMatchBegin; }      
  kUTRMS.kUnitigLength = atoi (flds[8]);				
  kUTRMS.readLength = atoi (flds[9]);				
  kUTRMS.kUnitigNumber = atoi (flds[10]); 
}

int processKUnitigVsReadMatches (char *readVsKUnitigFile, char* outputFileName)
{
     char *cptr;
     charb cmd(500), line(2000);
     FILE *infile;
     int numFlds;
   
     sprintf (cmd, "zcat -f %s", readVsKUnitigFile);
     infile = Popen (cmd, "r");
     if (! fgets (line, 2000, infile))
	  return (1); // A critical file doesn't exist
     outputFile=Fopen(outputFileName,"w");
     // Load the appropriate stuff
     numFlds = getFldsFromLine(line);
     cptr = flds[numFlds-1];
     rdPrefixHold[0] = cptr[0];
     rdPrefixHold[1] = cptr[1];
     cptr += 2;
     readNum = readNumHold = atoll (cptr);
     evenReadMatchStructs.clear();
     oddReadMatchStructs.clear();
     updateMatchRecords(readNum, cptr, flds);

     while (fgets (line, 2000, infile)) {
	  numFlds = getFldsFromLine(line);
	  cptr = flds[numFlds-1];
	  rdPrefix[0] = cptr[0];
	  rdPrefix[1] = cptr[1];
	  cptr += 2;
	  readNum = atof (cptr);
	  if ((strcmp (rdPrefix, rdPrefixHold) == 0) &&
	      readNum == readNumHold) {
	       // load more data
	       updateMatchRecords(readNum, cptr, flds);
	       continue;
	  }
	  if ((strcmp (rdPrefix, rdPrefixHold) != 0) ||
	      (readNum != readNumHold+1) ||
	      (readNum % 2 == 0)) {
	       // Get the super-read for the insert we just finished reading
	       getSuperReadsForInsert();
	       // Set up and load the new data
               evenReadMatchStructs.clear();
               oddReadMatchStructs.clear();
	       updateMatchRecords(readNum, cptr, flds);
	       // Update what the old data is
	       strcpy (rdPrefixHold, rdPrefix);
	       readNumHold = readNum;
	       continue;
	  }
	  // If we get here we've gotten to the second read of a mate pair
	  // load the data
	  updateMatchRecords(readNum, cptr, flds);
	  // hold the updated read info
	  readNumHold = readNum;
     }
     fclose (infile);
     // Output the stuff for the old pair
     getSuperReadsForInsert();
     fclose(outputFile);
     return (0);
}

// returns 1 if successful, 0 if too many nodes (so failure)
int joinKUnitigsFromMates (int insertLengthMean, int insertLengthStdev)
{
     
     int lastOffsetToTest = 6000, lastOffsetToTestIfNotMate2, maxOffsetToAllow;
     int j;
     struct unitigLocStruct unitigLocVal;
     struct abbrevUnitigLocStruct abbrevUnitigLocVal, *abbRLPtr;
     size_t maxNodes;
     char *vptr;
     int unitig1, unitig2;
     char ori; 
     int offset;
     int elementIndex;
     int overlapLength;
     int ahg, bhg;


     lastOffsetToTest = insertLengthMean+4*insertLengthStdev;
     // The following assumes that all the overlaps are of length
     // minOverlapLength
     lastOffsetToTestIfNotMate2 = lastOffsetToTest - (unitigLengths[mateUnitig2]-minOverlapLength);
     // Adjust overlaps for mateUnitig1 if ori not 'F' and for
     //    mateUnitig2 if ori not 'R' to what they would be if they had
     //    the desired orientation. Note that mateUnitig1 is only a
     //    source and mateUnitig2 is only a sink
     unitigLocVal.unitig2 = mateUnitig1;
     unitigLocVal.frontEdgeOffset = unitigLengths[mateUnitig1];
     unitigLocVal.ori = mateUnitig1ori;
     numTreesUsed = 0;
     abbrevUnitigLocVal.frontEdgeOffset = unitigLocVal.frontEdgeOffset;
     abbrevUnitigLocVal.ori = unitigLocVal.ori;
     // We may want to change the following
     if (treeArr[mateUnitig1].root == TREE_NIL)
     {
	  treeReinitList[numTreesUsed] = mateUnitig1;
	  numTreesUsed++;
     }
     abbrevUnitigLocVal.pathNum = 0;
     RBTreeInsertElement (treeArr + mateUnitig1, (char *) &abbrevUnitigLocVal);
#if 0
     printf ("Inserting in the RB tree at %d: fEO = %d, pN = %u ori = %c\n", mateUnitig1, abbrevUnitigLocVal.frontEdgeOffset, abbrevUnitigLocVal.pathNum, abbrevUnitigLocVal.ori);
#endif     
     forward_path_unitigs.push(unitigLocVal);
     
     startingNodes.clear();
     startingNodes.insert(unitigLocVal);
     nodeArray.clear();
     nodeToIndexMap.clear();
//     fwdEdgeList.clear();
//     revEdgeList.clear();
     nodeArray.push_back(unitigLocVal);
     nodeToIndexMap.insert (std::pair<unitigLocStruct, int> (unitigLocVal, nodeArray.size()-1) );
     maxNodes = 1;
//     printf ("Got to 30\n");
     while (!forward_path_unitigs.empty())
     {
       if (forward_path_unitigs.size() > maxNodes)
         maxNodes = forward_path_unitigs.size();
#if DEBUG
	  printf ("Starting new offset: "); fflush (stdout);
	  for(min_heap::iterator it = forward_path_unitigs.begin(); it != forward_path_unitigs.end(); ++it)
            printf("%d ", it->frontEdgeOffset); fflush(stdout);
	  printf ("\n"); fflush (stdout);
#endif
          unitigLocVal = forward_path_unitigs.pop();
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
	       treeFindElement (treeArr + unitig1, (char *) &abbrevUnitigLocVal);
	  assert (elementIndex != TREE_NIL);
	  setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	  abbRLPtr = (abbrevUnitigLocStruct *) vptr;
//	  printf ("Got to 40\n");
	  for (j = startOverlapByUnitig[unitig1];
	       j < startOverlapByUnitig[unitig1 + 1]; j++)
	  {
#if DEBUG
	       printf ("Starting an overlap:");
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
	       if (ori == 'F')
	       {
		    if (overlapData[j].bhg <= 0) 
			 continue;
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
			 continue;
		    ahg = overlapData[j].ahg;
		    if (overlapData[j].ori == 'N')
			 abbrevUnitigLocVal.ori = 'R';
		    else
			 abbrevUnitigLocVal.ori = 'F';
		    abbrevUnitigLocVal.frontEdgeOffset = offset - ahg;
	       }
	       // Skip if the offset is too large
#if DEBUG
	       printf ("frontEdgeOffset = %d, lastOffsetToTest = %d\n", abbrevUnitigLocVal.frontEdgeOffset, lastOffsetToTest);
#endif
	       if (unitig2 == mateUnitig2)
		    maxOffsetToAllow = lastOffsetToTest;
	       else
		    maxOffsetToAllow = lastOffsetToTestIfNotMate2;
	       if (abbrevUnitigLocVal.frontEdgeOffset > maxOffsetToAllow)
		    continue;
#if DEBUG
	       printf ("cur front = %d\n", abbrevUnitigLocVal.frontEdgeOffset);
#endif
	       // Skip if abbrevUnitigLocVal al unitig y seen for unitig2
	       elementIndex =
		    treeFindElement (treeArr + unitig2, (char *) &abbrevUnitigLocVal);
	       setTreeValPtr (vptr, treeArr + unitig2, elementIndex);
	       abbRLPtr = (abbrevUnitigLocStruct *) vptr;
	       
//	       printf ("Got to 60\n");
	       if (elementIndex != TREE_NIL)
		    continue;
#if DEBUG
	       printf ("Adding to the tree\n");
#endif
	       // Insert this value in the priority queue
	       unitigLocVal.unitig2 = unitig2;
	       unitigLocVal.frontEdgeOffset = abbrevUnitigLocVal.frontEdgeOffset;
	       unitigLocVal.ori = abbrevUnitigLocVal.ori;
               forward_path_unitigs.push(unitigLocVal);
	       if (treeArr[unitig2].root == TREE_NIL)
	       {
		    // If unitig2's tree never seen before
		    //    Add unitig2 to list of trees to reinit
		    treeReinitList[numTreesUsed] = unitig2;
		    numTreesUsed++;
	       }
	       // Add offset to list for tree
	       abbrevUnitigLocVal.pathNum = 0;
	       RBTreeInsertElement (treeArr + unitig2, (char *) &abbrevUnitigLocVal);
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
     }			// Ends !forward_path_unitigs.empty() line
     if (maxNodes > MAX_NODES_ALLOWED)
	  return (0);
     else
	  return (1);
}

// The following returns the overlap length if it is greater than the
// existing largest overlap on the end. It returns -1 if not.
int getOvlLenFromOvlIndicesPlus (int maxOvlIndex, int j, int maxOvlLen, int whichEnd)
{
     int thisUnitig, otherUnitig, overlapLength;
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
     if (ptr->ori == 'F') {
	  endOffset = ptr->frontEdgeOffset;
	  beginOffset = endOffset - unitigLengths[ptr->unitig1];
     }
     else {
	  beginOffset = ptr->frontEdgeOffset;
	  endOffset = beginOffset - unitigLengths[ptr->unitig1];
     }
//     fprintf (outfile, "uni = %d, offset = %d, ori = %c, beginOffset = %d, endOffset = %d, numOvlsIn = %d, numOvlsOut = %d\n", ptr->unitig1, ptr->frontEdgeOffset, ptr->ori, beginOffset, endOffset, ptr->numOverlapsIn, ptr->numOverlapsOut);
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].unitig1 = ptr->unitig1;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].frontEdgeOffset = ptr->frontEdgeOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].ori = ptr->ori;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].beginOffset = beginOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].endOffset = endOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = ptr->numOverlapsIn;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = ptr->numOverlapsOut;
     ++numUnitigPathPrintRecsOnPath;
     if (ptr->numOverlapsIn == 0)
	  approxNumPaths += ptr->numOverlapsOut;
     else if (ptr->numOverlapsOut > 0)
	  approxNumPaths += (ptr->numOverlapsOut - 1);
}

void completePathPrint (struct abbrevUnitigLocStruct *ptr)
{
     struct abbrevUnitigLocStruct abbrevUnitigLocVal, *abbRLPtr;
     struct unitigLocStruct unitigLocVal;
     struct unitigPathPrintStruct unitigPathPrintVal, *rppvPtr1, *rppvPtr2;
     struct unitigConnectionsForPathStruct unitigConnectionsForPathRec;
     char ori;
#ifdef NEW_STUFF
     char tempOri1, tempOri2;
#endif
     int isSpecialCase, finalOffset, minConnectingOffset, i, index;
     int unitig1, unitig2, offset, elementIndex1, elementIndex, overlapLength;
     int minConnectingUnitig=0, minConnectingOverlapIndex;
     char minConnectingOri=' ';
     char *vptr;
     double numStdevsFromMean;
     // In the following we assume we move from left to right when moving from
     // beginUnitig to endUnitig.
     ++curPathNum;
#if 0
     printf ("curPathNum = %d, nodePathNum = %d; ", curPathNum, ptr->pathNum);
     printf ("frontEdgeOffset = %d, ori = %c\n", ptr->frontEdgeOffset, ptr->ori);
#endif
     ori = ptr->ori;
     if (ori != endUnitigOri) return;
//     printf ("In rtn completePathPrint\n");
     isSpecialCase = 1;
     finalOffset = ptr->frontEdgeOffset;
     minConnectingOffset = finalOffset + 1000000;
     numStdevsFromMean = (finalOffset - insertLengthMeanBetweenKUnisForInsertGlobal)/insertLengthStdevGlobal;
#ifdef KILLED111115
     fprintf (outfile, "%d %f\n", finalOffset, numStdevsFromMean);
#endif
     for (i=startOverlapIndexByUnitig2[endUnitig]; i<startOverlapIndexByUnitig2[endUnitig+1]; i++) {
	  index = unitig2OverlapIndex[i];
	  unitig1 = overlapData[index].unitig1;
	  if (overlapData[index].ahg >= 0)
	       overlapLength = unitigLengths[unitig1] - overlapData[index].ahg;
	  else
	       overlapLength = unitigLengths[endUnitig] + overlapData[index].ahg;
	  if (overlapLength > unitigLengths[unitig1])
	       overlapLength = unitigLengths[unitig1];
	  if (overlapLength > unitigLengths[endUnitig])
	       overlapLength = unitigLengths[endUnitig];
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
	  elementIndex = treeFindElement (treeArr + unitig1, (char *) &abbrevUnitigLocVal);
	  if (elementIndex == TREE_NIL) continue;
//	  printf ("We are at 2\n");
	  setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	  abbRLPtr = (abbrevUnitigLocStruct *) vptr;
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
//	  fprintf (stderr, "We shouldn't get here\n");
	  fprintf (stdout, "We shouldn't get here\n");
	  unitigLocVal.unitig2 = minConnectingUnitig;
	  unitigLocVal.frontEdgeOffset = minConnectingOffset;
	  unitigLocVal.ori = minConnectingOri;
	  // ..and for the path
	  unitigPathPrintVal.unitig1 = endUnitig;
	  unitigPathPrintVal.numOverlapsIn = 1;
	  unitigPathPrintVal.numOverlapsOut = 0;
	  unitigPathPrintVal.ori = 'R'; // Forced; may be adjusted later
	  unitigPathPrintVal.frontEdgeOffset = finalOffset;
	  RBTreeInsertElement (treeArr2, (char *) &unitigPathPrintVal);
	  unitigPathPrintVal.unitig1 = minConnectingUnitig;
	  unitigPathPrintVal.frontEdgeOffset = minConnectingOffset;	
  unitigPathPrintVal.numOverlapsIn = 0;
	  unitigPathPrintVal.numOverlapsOut = 1;
	  unitigPathPrintVal.ori = minConnectingOri;
	  RBTreeInsertElement (treeArr2, (char *) &unitigPathPrintVal);
     }
     else {
	  unitigLocVal.unitig2 = endUnitig;
	  unitigLocVal.frontEdgeOffset = finalOffset;
	  unitigLocVal.ori = endUnitigOri;
	  unitigPathPrintVal.unitig1 = endUnitig;
	  unitigPathPrintVal.frontEdgeOffset = finalOffset;
	  unitigPathPrintVal.numOverlapsIn = 0;
	  unitigPathPrintVal.numOverlapsOut = 0;
	  unitigPathPrintVal.ori = endUnitigOri;
#if 0
	  printf ("Inserting unitig1 = %d, offset = %d, ori = %c\n", unitigPathPrintVal.unitig1, unitigPathPrintVal.frontEdgeOffset, unitigPathPrintVal.ori);
#endif
	  RBTreeInsertElement (treeArr2, (char *) &unitigPathPrintVal);
     }
#if 0
     printf ("isSpecialCase = %d, unitigLocVal = %d, %d, %c\n", isSpecialCase, endUnitig, finalOffset, unitigLocVal.ori);
#endif
     backward_path_unitigs.push(unitigLocVal);
     while (!backward_path_unitigs.empty()) {
          unitigLocVal = backward_path_unitigs.pop();
	  unitig2 = unitigLocVal.unitig2;
	  offset = unitigLocVal.frontEdgeOffset;
	  ori = unitigLocVal.ori;
	  unitigPathPrintVal.unitig1 = unitig2;
	  unitigPathPrintVal.frontEdgeOffset = offset;
	  unitigPathPrintVal.ori = ori;
	  elementIndex1 = treeFindElement (treeArr2, (char *) &unitigPathPrintVal);
	  setTreeValPtr (vptr, treeArr2, elementIndex1);
	  rppvPtr1 = (unitigPathPrintStruct *) vptr;
#if 0
	  printf ("unitig2 = %d, offset = %d, ori = %c; elementIndex1 = %d\n", unitig2, offset, ori, elementIndex1);
#endif
	  for (i=startOverlapIndexByUnitig2[unitig2]; i<startOverlapIndexByUnitig2[unitig2+1]; i++) {
	       index = unitig2OverlapIndex[i];
	       unitig1 = overlapData[index].unitig1;
#if 0
	       printf ("unitig1 = %d\n", unitig1);
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
	       elementIndex = treeFindElement (treeArr + unitig1, (char *) &abbrevUnitigLocVal);
//	       printf ("Got to 21, unitig1 = %d\n", unitig1);
//	       printf ("fEO = %d, pN = %u ori = %c\n", abbrevUnitigLocVal.frontEdgeOffset, abbrevUnitigLocVal.pathNum, abbrevUnitigLocVal.ori);
	       if (elementIndex == TREE_NIL) continue;
//	       printf ("Got to 215\n");
	       setTreeValPtr (vptr, treeArr + unitig1, elementIndex);
	       abbRLPtr = (abbrevUnitigLocStruct *) vptr;
//	       printf ("Got to 22\n");
	       if (abbRLPtr->frontEdgeOffset >= offset) continue;
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
                    backward_path_unitigs.push(unitigLocVal);
		    unitigPathPrintVal.unitig1 = unitig1;
		    unitigPathPrintVal.frontEdgeOffset = unitigLocVal.frontEdgeOffset;
		    unitigPathPrintVal.ori = unitigLocVal.ori;
		    unitigPathPrintVal.numOverlapsIn = 0;
		    unitigPathPrintVal.numOverlapsOut = 0;
		    RBTreeInsertElement (treeArr2, (char *) &unitigPathPrintVal);
	       }
		    
	       unitigPathPrintVal.unitig1 = unitig1;
	       unitigPathPrintVal.frontEdgeOffset = abbRLPtr->frontEdgeOffset;
	       unitigPathPrintVal.ori = abbRLPtr->ori;
	       elementIndex = treeFindElement (treeArr2, (char *) &unitigPathPrintVal);
	       setTreeValPtr (vptr, treeArr2, elementIndex);
	       rppvPtr2 = (unitigPathPrintStruct *) vptr;
	       int frontEdgeOffset1 = rppvPtr2->frontEdgeOffset;
	       // The following must be recalced in case the array had to be
	       // moved due to needing more space
	       setTreeValPtr (vptr, treeArr2, elementIndex1);
	       rppvPtr1 = (unitigPathPrintStruct *) vptr;
	       int frontEdgeOffset2 = rppvPtr1->frontEdgeOffset;
	       if (frontEdgeOffset2 - frontEdgeOffset1 != unitigLengths[rppvPtr1->unitig1] - minOverlapLength)
		    continue;
	       setTreeValPtr (vptr, treeArr2, elementIndex);
	       rppvPtr2 = (unitigPathPrintStruct *) vptr;
	       ++(rppvPtr2->numOverlapsOut);
	       if (rppvPtr2->numOverlapsOut > 1) {
		    if (rppvPtr2->frontEdgeOffset < splitJoinWindowMin)
			 splitJoinWindowMin = rppvPtr2->frontEdgeOffset; }
	       // The following must be recalced in case the array had to be
	       // moved due to needing more space
	       setTreeValPtr (vptr, treeArr2, elementIndex1);
	       rppvPtr1 = (unitigPathPrintStruct *) vptr;
	       ++(rppvPtr1->numOverlapsIn);
	       if (rppvPtr1->numOverlapsIn > 1) {
		    if (rppvPtr1->frontEdgeOffset > splitJoinWindowMax)
			 splitJoinWindowMax = rppvPtr1->frontEdgeOffset; }
#ifdef NEW_STUFF
	       tempOri1 = rppvPtr1->ori;
	       if (rppvPtr1->unitig1 == beginUnitig) tempOri1 = beginUnitigOri;
	       if (rppvPtr1->unitig1 == endUnitig) tempOri1 = endUnitigOri;
	       tempOri2 = rppvPtr2->ori;
	       if (rppvPtr2->unitig1 == beginUnitig) tempOri2 = beginUnitigOri;
	       if (rppvPtr2->unitig1 == endUnitig) tempOri2 = endUnitigOri;
//	       printf ("Node (%d, %d, %c) -> (%d, %d, %c)\n", rppvPtr2->unitig1, rppvPtr2->frontEdgeOffset, tempOri2, rppvPtr1->unitig1, rppvPtr1->frontEdgeOffset, tempOri1);
	       unitigConnectionsForPathRec.unitig1 = rppvPtr2->unitig1;
	       unitigConnectionsForPathRec.unitig2 = rppvPtr1->unitig1;
	       unitigConnectionsForPathRec.frontEdgeOffset1 = rppvPtr2->frontEdgeOffset;
	       unitigConnectionsForPathRec.frontEdgeOffset2 = rppvPtr1->frontEdgeOffset;
	       unitigConnectionsForPathRec.ori1 = tempOri2;
	       unitigConnectionsForPathRec.ori2 = tempOri1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].unitig1 = unitigConnectionsForPathRec.unitig1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].unitig2 = unitigConnectionsForPathRec.unitig2;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].frontEdgeOffset1 = unitigConnectionsForPathRec.frontEdgeOffset1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].frontEdgeOffset2 = unitigConnectionsForPathRec.frontEdgeOffset2;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].ori1 = unitigConnectionsForPathRec.ori1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].ori2 = unitigConnectionsForPathRec.ori2;
//	       unitigConnectionsForPathData.push_back (unitigConnectionsForPathRec);
	       ++numUnitigConnectionsForPathData;
#endif
//	       printf ("Got to 24\n");
	       
	  }
     }
// #ifdef KILLED111115
     for (i=0; i<(int) unitigConnectionsForPathData.size(); i++) {
	  unitigLocStruct uLS;
//	  int index1, index2;
	  int firstIndex, secondIndex;
	  std::map<unitigLocStruct, int>::iterator it;
	  std::list<int>::iterator l_it;
	  struct nodePair nodeValue;
	  uLS.unitig2 = unitigConnectionsForPathData[i].unitig1;
	  uLS.frontEdgeOffset = unitigConnectionsForPathData[i].frontEdgeOffset1;
	  uLS.ori = unitigConnectionsForPathData[i].ori1;
	  it = nodeToIndexMap.find (uLS);
	  if (it == nodeToIndexMap.end()) {
	       nodeArray.push_back(uLS);
	       firstIndex = nodeArray.size()-1;
	       nodeToIndexMap.insert(std::pair<unitigLocStruct, int> (uLS, firstIndex ) ); }
	  else {
	       firstIndex = it->second; }
	  uLS.unitig2 = unitigConnectionsForPathData[i].unitig2;
	  uLS.frontEdgeOffset = unitigConnectionsForPathData[i].frontEdgeOffset2;
	  uLS.ori = unitigConnectionsForPathData[i].ori2;
	  it = nodeToIndexMap.find (uLS);
	  if (it == nodeToIndexMap.end()) {
	       nodeArray.push_back(uLS);
	       secondIndex = nodeArray.size()-1;
	       nodeToIndexMap.insert(std::pair<unitigLocStruct, int> (uLS, secondIndex ) ); }
	  else {
	       secondIndex = it->second; }
	  nodeValue.node1 = firstIndex;
	  nodeValue.node2 = secondIndex;
	  edgeList.insert(std::pair<int, int> (firstIndex, secondIndex));
	  
#ifdef KILL120102	  
	  fprintf (stderr, "%s%lld Node (%d, %d, %c) -> (%d, %d, %c)\n", rdPrefixHold, readNumHold, unitigConnectionsForPathData[i].unitig1, unitigConnectionsForPathData[i].frontEdgeOffset1, unitigConnectionsForPathData[i].ori1, unitigConnectionsForPathData[i].unitig2, unitigConnectionsForPathData[i].frontEdgeOffset2, unitigConnectionsForPathData[i].ori2);
#endif
     }
// #endif
#if 0
     printf ("tree root = %d\n", treeArr2[0].root);
#endif
     numUnitigPathPrintRecsOnPath = 0;
     if (treeSize <= maxDiffInsertSizesForPrinting)
	  inOrderTreeWalk (treeArr2, treeArr2[0].root,
			   (void (*)(char *)) printPathNode);
#ifdef KILLED111115
     for (i=0; i<numUnitigPathPrintRecsOnPath; i++)
	  fprintf (outfile, "uni = %d, offset = %d, ori = %c, beginOffset = %d, endOffset = %d, numOvlsIn = %d, numOvlsOut = %d\n", augmentedUnitigPathPrintData[i].unitig1, augmentedUnitigPathPrintData[i].frontEdgeOffset, augmentedUnitigPathPrintData[i].ori, augmentedUnitigPathPrintData[i].beginOffset, augmentedUnitigPathPrintData[i].endOffset, augmentedUnitigPathPrintData[i].numOverlapsIn, augmentedUnitigPathPrintData[i].numOverlapsOut);
#endif
     if (approxNumPaths == 1) {
	  generateSuperReadPlacementLinesForJoinedMates();
#if 0
	  int isReversed, superReadLength;
	  charb tempOutputString(100);
	  // the following uses augmentedUnitigPathPrintData
	  superReadLength = getSuperReadLength ();
	  isReversed = setSuperReadNameFromAugmentedPath ();
	  sprintf (outputString,"%s%lld %s ", rdPrefixHold, readNumHold-1, superReadName);
	  if (! isReversed)
	       sprintf (tempOutputString,"%d F\n", lengthAdjustment1);
	  else
	       sprintf (tempOutputString,"%d R\n", superReadLength - lengthAdjustment1);
	  strcat (outputString, tempOutputString);
	  sprintf (tempOutputString,"%s%lld %s ", rdPrefixHold, readNumHold, superReadName);
	  strcat (outputString, tempOutputString);
	  if (! isReversed)
	       sprintf (tempOutputString,"%d R\n", superReadLength - lengthAdjustment2);
	  else
	       sprintf (tempOutputString,"%d F\n", lengthAdjustment2);
	  strcat (outputString, tempOutputString);
#endif
     }
#if 0
     printf ("final offset = %d, arraySize = %d\n", finalOffset, dataArr2.arraySize);
#endif
     treeArr2[0].root = TREE_NIL;
     dataArr2.arraySize = 0;
}

void generateSuperReadPlacementLinesForJoinedMates (void)
{
     int isReversed, superReadLength;
     charb tempOutputString(100);
     // the following uses augmentedUnitigPathPrintData
     superReadLength = getSuperReadLength ();
     isReversed = setSuperReadNameFromAugmentedPath ();
     sprintf (outputString,"%s%lld %s ", rdPrefixHold, readNumHold-1, superReadName);
     if (! isReversed)
	  sprintf (tempOutputString,"%d F\n", lengthAdjustment1);
     else
	  sprintf (tempOutputString,"%d R\n", superReadLength - lengthAdjustment1);
     strcat (outputString, tempOutputString);
     sprintf (tempOutputString,"%s%lld %s ", rdPrefixHold, readNumHold, superReadName);
     strcat (outputString, tempOutputString);
     if (! isReversed)
	  sprintf (tempOutputString,"%d R\n", superReadLength - lengthAdjustment2);
     else
	  sprintf (tempOutputString,"%d F\n", lengthAdjustment2);
     strcat (outputString, tempOutputString);
}

int setSuperReadNameFromAugmentedPath (void)
{
     int isReversed=0, i;
     char *cptr;
     for (i=0; i<numUnitigPathPrintRecsOnPath/2; i++) {
	  if (augmentedUnitigPathPrintData[i].unitig1 != augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].unitig1) {
	       if (augmentedUnitigPathPrintData[i].unitig1 < augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].unitig1)
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (augmentedUnitigPathPrintData[i].ori == augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].ori) {
	       if (augmentedUnitigPathPrintData[i].ori == 'F')
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
     }
     cptr = superReadName;
     if (isReversed == 0) {
	  sprintf (cptr, "%d%c", augmentedUnitigPathPrintData[0].unitig1, augmentedUnitigPathPrintData[0].ori);
	  cptr += strlen (cptr);
	  for (i=1; i<numUnitigPathPrintRecsOnPath; i++) {
//	       sprintf (cptr, "_%d_%d%c", minOverlapLength, augmentedUnitigPathPrintData[i].unitig1, augmentedUnitigPathPrintData[i].ori);
	       sprintf (cptr, "_%d%c", augmentedUnitigPathPrintData[i].unitig1, augmentedUnitigPathPrintData[i].ori);
	       cptr += strlen (cptr);
	  }
     }
     else {
	  sprintf (cptr, "%d%c", augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-1].unitig1, (augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-1].ori == 'F') ? 'R' : 'F');
	  cptr += strlen (cptr);
	  for (i=numUnitigPathPrintRecsOnPath-2; i>=0; i--) {
//	       sprintf (cptr, "_%d_%d%c", minOverlapLength, augmentedUnitigPathPrintData[i].unitig1, (augmentedUnitigPathPrintData[i].ori == 'F') ? 'R' : 'F');
	       sprintf (cptr, "_%d%c", augmentedUnitigPathPrintData[i].unitig1, (augmentedUnitigPathPrintData[i].ori == 'F') ? 'R' : 'F');
	       cptr += strlen (cptr);
	  }
     }
     return (isReversed);
}

int getSuperReadLength(void)
{
     int totLen, i;
     totLen = unitigLengths[augmentedUnitigPathPrintData[0].unitig1];
     for (i=1; i<numUnitigPathPrintRecsOnPath; i++)
	  totLen += (unitigLengths[augmentedUnitigPathPrintData[i].unitig1] - minOverlapLength);
     
     return (totLen);
}

void funcToGetTreeSize (abbrevUnitigLocStruct *ptr)
{
     struct unitigLocStruct localUnitigLoc;
     if (ptr->ori == endUnitigOri) {
	  ++treeSize;
	  localUnitigLoc.unitig2 = endUnitig;
	  localUnitigLoc.frontEdgeOffset = ptr->frontEdgeOffset;
	  localUnitigLoc.ori = endUnitigOri;
	  endingNodes.insert (localUnitigLoc); 
	  nodeArray.push_back(localUnitigLoc);
	  nodeToIndexMap.insert (std::pair<unitigLocStruct, int> (localUnitigLoc, nodeArray.size()-1) ); }
}

bool unitigLocStructCompare (struct unitigLocStruct uLS1,
			    struct unitigLocStruct uLS2)
{
     if (uLS1.frontEdgeOffset != uLS2.frontEdgeOffset)
	  return (uLS1.frontEdgeOffset < uLS2.frontEdgeOffset);
     if (uLS1.unitig2 != uLS2.unitig2)
	  return (uLS1.unitig2 < uLS2.unitig2);
     return (uLS1.ori < uLS2.ori);
}

bool unitigLocStructCompareReversed (struct unitigLocStruct uLS1,
				    struct unitigLocStruct uLS2)
{
     if (uLS1.frontEdgeOffset < uLS2.frontEdgeOffset)
	  return (-1);
     if (uLS1.frontEdgeOffset > uLS2.frontEdgeOffset)
	  return (1);
     if (uLS1.unitig2 < uLS2.unitig2)
	  return (-1);
     if (uLS1.unitig2 > uLS2.unitig2)
	  return (1);
     if (uLS1.ori < uLS2.ori)
	  return (-1);
     if (uLS1.ori > uLS2.ori)
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

int getInt (const char *fname)
{
     FILE *infile;
     int tval;

     infile = Fopen (fname, "r");
     if (! fscanf (infile, "%d\n", &tval)) {
	  fprintf (stderr, "Couldn't read file %s. Bye!\n", fname);
	  exit (1);
     }
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

void findSingleReadSuperReads(char *readName)
{
     long long tempInt;
     char *cptr=readName+2;
     int countOfMatchingKUnitigs, offsetOfReadInSuperRead;
     int minReadOffset, maxReadOffset, minReadOffsetSeen, maxReadOffsetSeen;
     int i, j, recNumToUse=0;
     int isReversed=0;
     struct kuniToReadMatchStruct *kUTRMSptr;

     tempInt = atoll(cptr);
#ifdef KILLED111115
     printf ("findSingleReadSuperReads\n");
#endif
     if (tempInt % 2 == 0) {
          countOfMatchingKUnitigs = evenReadMatchStructs.size();
	  kUTRMSptr = &(evenReadMatchStructs[0]);
     }
     else {
          countOfMatchingKUnitigs = oddReadMatchStructs.size();
	  kUTRMSptr = &(oddReadMatchStructs[0]);
     }
     
//     printf ("countOfMatchingKUnitigs = %d\n", countOfMatchingKUnitigs);
     i = 0;
     minReadOffsetSeen = kUTRMSptr[i].readMatchBegin;
     maxReadOffsetSeen = kUTRMSptr[i].readMatchEnd;
     matchStructIsUsed[i] = 1;
     for (i=1; i<countOfMatchingKUnitigs; i++) {
	  matchStructIsUsed[i] = 0;
	  if (kUTRMSptr[i].readMatchEnd <= maxReadOffsetSeen)
	       continue;
	  if (kUTRMSptr[i].readMatchBegin < maxReadOffsetSeen-minOverlapLength)
	       continue; // Otherwise the k-unitigs overlap too much
	  if (kUTRMSptr[i].readMatchBegin > maxReadOffsetSeen)
	       return; // Part of the middle of the read is uncovered by k-unis
	  maxReadOffsetSeen = kUTRMSptr[i].readMatchEnd;
	  matchStructIsUsed[i] = 1;
//	  printf ("Struct number %d is used, kUnitig = %d\n", i, kUTRMSptr[i].kUnitigNumber);
     }
     if (minReadOffsetSeen + (kUTRMSptr[0].readLength - maxReadOffsetSeen) > maxTotAllowableMissingOnEnds)
	  return;
     
     i=-1; j=countOfMatchingKUnitigs;
     isReversed = 0;
     while (1) {
	  ++i; --j;
	  while (!matchStructIsUsed[i])
	       ++i;
	  while (!matchStructIsUsed[j])
	       --j;
	  if (kUTRMSptr[i].kUnitigNumber != kUTRMSptr[j].kUnitigNumber) {
	       if (kUTRMSptr[i].kUnitigNumber < kUTRMSptr[j].kUnitigNumber)
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (kUTRMSptr[i].ori == kUTRMSptr[j].ori) {
	       if (kUTRMSptr[i].ori == 'F')
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (j<=i)
	       break;
     }
     cptr = superReadName;
     if (isReversed == 0) {
	  for (i=0; 1; i++)
	       if (matchStructIsUsed[i])
		    break;
	  recNumToUse = i;
	  sprintf (cptr, "%d%c", kUTRMSptr[i].kUnitigNumber, kUTRMSptr[i].ori);
	  cptr += strlen(cptr);
	  maxReadOffset = kUTRMSptr[i].readMatchEnd;
	  for (++i; i<countOfMatchingKUnitigs; i++) {
	       if (! matchStructIsUsed[i])
		    continue;
	       // The next is the overlap amount between k-unitigs, which we now require to be minOverlapLength
	       if (maxReadOffset-kUTRMSptr[i].readMatchBegin != minOverlapLength)
		    return;
//	       sprintf (cptr, "_%d_%d%c", maxReadOffset-kUTRMSptr[i].readMatchBegin, kUTRMSptr[i].kUnitigNumber, kUTRMSptr[i].ori);
	       sprintf (cptr, "_%d%c", kUTRMSptr[i].kUnitigNumber, kUTRMSptr[i].ori);
	       cptr += strlen(cptr);
	       maxReadOffset = kUTRMSptr[i].readMatchEnd;
	  }
	  // Must do the output here
	  if (kUTRMSptr[recNumToUse].ori == 'F')
	       offsetOfReadInSuperRead = kUTRMSptr[recNumToUse].kUnitigMatchBegin - kUTRMSptr[recNumToUse].readMatchBegin;
	  else
	       offsetOfReadInSuperRead = unitigLengths[kUTRMSptr[recNumToUse].kUnitigNumber] - kUTRMSptr[recNumToUse].kUnitigMatchEnd - kUTRMSptr[i].readMatchBegin;
//	  printf ("%s %s %d %c\n", readName, superReadName, offsetOfReadInSuperRead, kUTRMSptr[recNumToUse].ori);
	  fprintf (outputFile, "%s %s %d %c\n", readName, superReadName, offsetOfReadInSuperRead, 'F');	  
     }
     else { // The k-unitigs are reversed from those reported
	  for (i=countOfMatchingKUnitigs-1; 1; i--)
	       if (matchStructIsUsed[i])
		    break;
	  recNumToUse = i;
	  sprintf (cptr, "%d%c", kUTRMSptr[i].kUnitigNumber, (kUTRMSptr[i].ori == 'F') ? 'R' : 'F');
	  cptr += strlen(cptr);
	  minReadOffset = kUTRMSptr[i].readMatchBegin;
	  for (--i; i>=0; i--) {
	       if (! matchStructIsUsed[i])
		    continue;
	       // The next is the overlap amount between k-unitigs, which we now require to be minOverlapLength
	       if (kUTRMSptr[i].readMatchEnd-minReadOffset != minOverlapLength)
		    return;
//	       sprintf (cptr, "_%d_%d%c", kUTRMSptr[i].readMatchEnd-minReadOffset, kUTRMSptr[i].kUnitigNumber, (kUTRMSptr[i].ori == 'F') ? 'R' : 'F');
	       sprintf (cptr, "_%d%c", kUTRMSptr[i].kUnitigNumber, (kUTRMSptr[i].ori == 'F') ? 'R' : 'F');
	       cptr += strlen(cptr);
	       minReadOffset = kUTRMSptr[i].readMatchBegin;
	  }
	  // Must do the output here
	  if (kUTRMSptr[recNumToUse].ori == 'F')
	       offsetOfReadInSuperRead = (unitigLengths[kUTRMSptr[recNumToUse].kUnitigNumber]-kUTRMSptr[recNumToUse].kUnitigMatchBegin) + kUTRMSptr[recNumToUse].readMatchBegin;
	  else
	       offsetOfReadInSuperRead = kUTRMSptr[recNumToUse].kUnitigMatchEnd + kUTRMSptr[recNumToUse].readMatchBegin;
	  // The k-unitigs are reversed from those reported
//	  printf ("%s %s %d %c\n", readName, superReadName, offsetOfReadInSuperRead, (kUTRMSptr[recNumToUse].ori == 'F') ? 'R' : 'F');
	  fprintf (outputFile, "%s %s %d %c\n", readName, superReadName, offsetOfReadInSuperRead, 'R');
     }
//     printf ("At 50\n");
}

void getSuperReadsForInsert (void)
{
     char readNameSpace[200];
     charb tempStderrStr(200);
     int insertLengthMean;
     int successCode;
     struct abbrevUnitigLocStruct abbULS1;
     struct unitigLocStruct tempULS;
     int elementIndex=0, elementIndex2=0;
     int numPossibleLengths=0;
     int startValue;
     int pathNum=0;
     std::set<unitigLocStruct>::iterator it1;
     std::map<unitigLocStruct, int>::iterator it2;
     int numNodes = 0;
     ExpandingBuffer<int> pathNumArray;
     std::stack<int> nodeIntArray;
     int localNodeNumber = 0, overlapMatchIndexHold = 0;
     int localUnitigNumber = 0, localNodeNumberHold = 0;
     int lastGoodNodeNumber;
     int localFrontEdgeOffset = 0, localSuperReadLength = 0;
     int doMinimalWorkHere;
     int distFromEndOfSuperRead = 0;

     // Output the stuff for the old pair
     stderrOutputString[0] = 0;
#ifdef KILLED111115
     printf ("%s%lld %ld %ld\n", rdPrefixHold, readNumHold, evenReadMatchStructs.size(), oddReadMatchStructs.size());
#endif
     sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
     if (evenReadMatchStructs.empty() || oddReadMatchStructs.empty()) {
	  findSingleReadSuperReads(readNameSpace);
	  return; }
     // If we get here both the even read and the odd read have
     // matches to k-unitigs
     // The next takes care of the case where both the source and
     // destination k-unitig are the same (We don't join in this case)
     mateUnitig1 = evenReadMatchStructs[0].kUnitigNumber;
     mateUnitig2 = oddReadMatchStructs[0].kUnitigNumber;
     if (mateUnitig1 == mateUnitig2) {
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold-1);
	  findSingleReadSuperReads (readNameSpace);
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
	  findSingleReadSuperReads (readNameSpace);
	  ++numPairsInOneUnitig;
	  return;
     }
     if ((evenReadMatchStructs[0].readMatchBegin + oddReadMatchStructs[0].readMatchBegin <= maxTotAllowableMissingOnEnds)) {
	  mateUnitig1ori = evenReadMatchStructs[0].ori;
	  if (oddReadMatchStructs[0].ori == 'F')
	       mateUnitig2ori = 'R';
	  else
	       mateUnitig2ori = 'F';
	  if (mateUnitig1ori == 'F')
	       lengthAdjustment1 = evenReadMatchStructs[0].ahg;
	  else
	       lengthAdjustment1 = - evenReadMatchStructs[0].bhg;
	  if (mateUnitig2ori == 'R')
	       lengthAdjustment2 = oddReadMatchStructs[0].ahg;
	  else
	       lengthAdjustment2 = - oddReadMatchStructs[0].bhg;
	  insertLengthMean = mean[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]] + (lengthAdjustment1 + lengthAdjustment2);
//	  printf ("lA1 = %d; lA2 = %d\n", lengthAdjustment1, lengthAdjustment2);
		   
#ifdef KILLED111115
	  printf ("joinKUnitigsFromMates for pair %s%lld %s%lld using mean %d\n", rdPrefixHold, readNumHold-1, rdPrefixHold, readNumHold, insertLengthMean);
#endif
	  // The following to check what we're doing
	  insertLengthMeanBetweenKUnisForInsertGlobal = insertLengthMean;
	  insertLengthStdevGlobal = stdev[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]];
	  successCode = joinKUnitigsFromMates (insertLengthMean, stdev[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]]);
	  if (! successCode)
	       goto afterSuperRead;
	  // Now doing the back trace
	  curPathNum = 0;
	  approxNumPaths = 0;
	  beginUnitig = mateUnitig1; beginUnitigOri = mateUnitig1ori;
	  endUnitig = mateUnitig2; endUnitigOri = mateUnitig2ori;
	  if (treeArr[mateUnitig2].root == TREE_NIL)
	       goto afterSuperRead;
	  treeSize = 0;
	  edgeList.clear();
	  endingNodes.clear();
	  fwdConnections.clear();
	  revConnections.clear();
	  inOrderTreeWalk (treeArr + endUnitig, treeArr[endUnitig].root,
			   (void (*)(char *)) funcToGetTreeSize);
#ifdef KILLED111115
	  printf ("treeSize = %d\n", treeSize);
#ifndef NO_OUTPUT
	  fprintf (outfile, "endUnitig = %d\n", endUnitig); // This prints
#endif
#endif
	  // This is where the main print statement is
	  splitJoinWindowMin = INT_MAX;
	  splitJoinWindowMax = INT_MIN;
	  // Don't bother if treeSize > 1; we must disambiguate
	  if (treeSize > 1) {
	       splitJoinWindowMin = INT_MIN;
	       splitJoinWindowMax = INT_MAX; }

	  unitigConnectionsForPathData.clear();
	  numUnitigConnectionsForPathData = 0;
	  inOrderTreeWalk (treeArr + endUnitig, treeArr[endUnitig].root,
			   (void (*)(char *)) completePathPrint);
	  if (approxNumPaths == 1)
	       ++numSimplyJoinable;
	  if (approxNumPaths <= 1)
	       goto afterSuperRead;
#if 1
	  sort (nodeArray.begin(), nodeArray.end(), unitigLocStructCompare);
	  for (unsigned int i=0; i<nodeArray.size(); i++) {
	       std::map<unitigLocStruct, int>::iterator it;
	       it = nodeToIndexMap.find (nodeArray[i]);
	       int j = it->second;
	       it->second = i;
	       newNodeNumsFromOld[j] = i; }
	  sortedEdgeList.clear();
	  for (std::set<std::pair<int, int> >::iterator it3141=edgeList.begin(); it3141 != edgeList.end(); it3141++)
	       sortedEdgeList.insert(std::pair<int, int> (newNodeNumsFromOld[it3141->first], newNodeNumsFromOld[it3141->second]));
	  edgeList = sortedEdgeList;
#ifdef KILL120102
	  fprintf (stderr, "Sorted node list:\n");
	  for (unsigned int i=0; i<nodeArray.size(); i++)
	       fprintf (stderr, "uni = %d, frontEdgeOffset = %d, ori = %c\n", (int) nodeArray[i].unitig2, (int) nodeArray[i].frontEdgeOffset, nodeArray[i].ori);
#endif
	  sprintf (stderrOutputString, "%s%lld\n", rdPrefixHold, readNumHold);
	  for (std::set<std::pair<int, int> >::iterator it3141=edgeList.begin(); it3141 != edgeList.end(); it3141++) {
	       sprintf (tempStderrStr, "Node %d %d %c -> %d %d %c\n", (int) nodeArray[it3141->first].unitig2, (int) nodeArray[it3141->first].frontEdgeOffset, (char) nodeArray[it3141->first].ori, (int) nodeArray[it3141->second].unitig2, (int) nodeArray[it3141->second].frontEdgeOffset, (char) nodeArray[it3141->second].ori);
	       strcat(stderrOutputString, tempStderrStr); }
#endif
	  struct nodePair tNodePair;
	  for (std::set<std::pair<int, int> >::iterator it3141=edgeList.begin(); it3141 != edgeList.end(); it3141++) {
	       tNodePair.node1 = it3141->first;
	       tNodePair.node2 = it3141->second;
	       fwdConnections.push_back(tNodePair);
	       revConnections.push_back(tNodePair);
	  }
	  sort (fwdConnections.begin(), fwdConnections.end(), firstNodeSort);
	  sort (revConnections.begin(), revConnections.end(), secondNodeSort);
	  numNodes = nodeArray.size();
	  fwdStartIndices.clear();
	  revStartIndices.clear();
	  fwdNumIndices.clear();
	  revNumIndices.clear();
	  pathNumArray.clear();
	  for (unsigned int j=0; j<nodeArray.size(); j++) {
	       fwdNumIndices[j] = revNumIndices[j] = 0;
	       fwdStartIndices[j] = revStartIndices[j] = 0; }
	  ++pathNum;
	  for (unsigned int j=0; j<nodeArray.size(); j++)
	       pathNumArray[j] = pathNum;
	  for (unsigned int j=0; j<fwdConnections.size(); j++) {
	       fwdStartIndices[fwdConnections[j].node1] = j;
	       ++fwdNumIndices[fwdConnections[j].node1];
	  }
#ifdef KILL120102
	  fprintf (stderr, "At 3 revStartIndices[1] = %d\n", (int) revStartIndices[1]);
#endif
	  for (unsigned int j=0; j<revConnections.size(); j++) {
	       revStartIndices[revConnections[j].node2] = j;
	       ++revNumIndices[revConnections[j].node2]; }
#ifdef KILL120102
	  for (unsigned int j=0; j<nodeArray.size(); j++)
	       fprintf (stderr, "At 4 revStartIndices[1] = %d\n", (int) revStartIndices[1]);
#endif
	  for (unsigned int j=0; j<fwdConnections.size(); j++) {
	       if (fwdNumIndices[j] > 0)
		    fwdStartIndices[j] -= (fwdNumIndices[j]-1);
	       else
		    fwdStartIndices[j] = 0; }
#ifdef KILL120102
	  fprintf (stderr, "At 5 revStartIndices[1] = %d\n", (int) revStartIndices[1]);
#endif
	  for (unsigned int j=0; j<revConnections.size(); j++) {
	       if (revNumIndices[j] > 0)
		    revStartIndices[j] -= (revNumIndices[j]-1);
	       else
		    revStartIndices[j] = 0; }
#ifdef KILL120102
	  fprintf (stderr, "revStartIndices[1] = %d\n", (int) revStartIndices[1]);
#endif
//     mustSplit1:
	  doMinimalWorkHere = 0;
	  if (evenReadMatchStructs.size() == 1) {
	       doMinimalWorkHere = 1; }
	  if (evenReadMatchStructs[0].ori == 'F')
	       startValue = evenReadMatchStructs[0].ahg;
	  else
	       startValue = - evenReadMatchStructs[0].bhg;
	  startValue += evenReadMatchStructs[0].readLength;
#ifdef KILL120102
	  fprintf (stderr, "Starting nodes:\n");
	  for (it1=startingNodes.begin(); it1!= startingNodes.end(); it1++)
	       fprintf (stderr, "%d %d %c\n", it1->unitig2, it1->frontEdgeOffset, it1->ori);
	  fprintf (stderr, "Ending nodes:\n");
	  for (it1=endingNodes.begin(); it1!= endingNodes.end(); it1++)
	       fprintf (stderr, "%d %d %c\n", it1->unitig2, it1->frontEdgeOffset, it1->ori);
#endif
	  for (int i=evenReadMatchStructs.size()-1; i>=0; i--) {
	       abbULS1.ori = evenReadMatchStructs[i].ori;
	       if (evenReadMatchStructs[i].ori == 'F')
		    abbULS1.frontEdgeOffset = startValue - evenReadMatchStructs[i].bhg;
	       else
		    abbULS1.frontEdgeOffset = startValue + evenReadMatchStructs[i].ahg;
	       elementIndex = treeFindElement (treeArr + evenReadMatchStructs[i].kUnitigNumber, (char *) &abbULS1);
	       if (elementIndex == TREE_NIL) {
#ifdef KILL120102
		    fprintf (stderr, "%s %d %d %d FAIL %d %d\n", readNameSpace, treeSize, i, abbULS1.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax);
#endif
		    continue; }
	       // If we get here we have a unitig on the path
	       char *vptr;
	       setTreeValPtr (vptr, treeArr+evenReadMatchStructs[i].kUnitigNumber, elementIndex);
	       pathNum = ((abbrevUnitigLocStruct *) vptr)->pathNum;
#ifdef KILL120102
	       fprintf (stderr, "%s %d %d %d SUCCESS %d %d %d\n", readNameSpace, treeSize, i, abbULS1.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax, pathNum);
#endif
	       if (pathNum == 0)
		    continue;
	       if (abbULS1.frontEdgeOffset >= splitJoinWindowMax)
		    continue;
	       // Checking for useless result
	       if (abbULS1.frontEdgeOffset <= splitJoinWindowMin) {
		    doMinimalWorkHere = 1;
		    break; }
	       // If we get here we've passed all the tests and we can use it
	       localUnitigNumber = evenReadMatchStructs[i].kUnitigNumber;
	       overlapMatchIndexHold = i;
	       break;
	  }
	  if (elementIndex == TREE_NIL) {
	       fprintf (stderr, "We should never get to TREE_NIL 1\n");
	       goto mustSplit2; }
	  if (doMinimalWorkHere) {
	       localUnitigNumber = evenReadMatchStructs[0].kUnitigNumber;
	       abbULS1.frontEdgeOffset = unitigLengths[localUnitigNumber];
	       abbULS1.ori = evenReadMatchStructs[0].ori; }
	  tempULS.unitig2 = localUnitigNumber;
	  tempULS.frontEdgeOffset = abbULS1.frontEdgeOffset;
	  tempULS.ori = abbULS1.ori;
#ifdef KILL120102
	  fprintf (stderr, "unitig2 = %d, frontEdgeOffset = %d, ori = %c\n", tempULS.unitig2, tempULS.frontEdgeOffset, tempULS.ori);
#endif
	  localNodeNumber = localNodeNumberHold = nodeToIndexMap[tempULS];
#ifdef KILL120102
	  fprintf (stderr, "localNodeNumber = %d\n", localNodeNumber);
#endif
	  ++pathNum;
	  pathNumArray[localNodeNumber] = pathNum;
	  nodeIntArray.push(localNodeNumber);
//	  fprintf (stderr, "fwdConnections\n");
//	  for (unsigned int j=0; j<fwdConnections.size(); j++)
//	       fprintf (stderr, "%d %d\n", fwdConnections[j].node1, fwdConnections[j].node2);
	  while (! nodeIntArray.empty()) {
	       int localLoopNodeNumber = nodeIntArray.top();
	       nodeIntArray.pop();
//	       fprintf (stderr, "localLoopNodeNumber = %d, start index = %d, num indices = %d\n", (int) localLoopNodeNumber, (int) fwdStartIndices[localLoopNodeNumber], (int) fwdNumIndices[localLoopNodeNumber]);
	       for (int j=fwdStartIndices[localLoopNodeNumber]; j<fwdStartIndices[localLoopNodeNumber]+fwdNumIndices[localLoopNodeNumber]; j++) {
		    localNodeNumber = fwdConnections[j].node2;
//		    fprintf (stderr, "pathNum = %d, localPathVal = %d\n", (int) pathNum, (int) pathNumArray[localNodeNumber]);
		    if (pathNumArray[localNodeNumber] < pathNum) {
			 nodeIntArray.push (localNodeNumber);
//			 fprintf (stderr, "localLoopNodeNumber = %d, fwdStartIndices = %d, numFwdIndices = %d, Pushing %d\n", localLoopNodeNumber, fwdStartIndices[localLoopNodeNumber], fwdNumIndices[localLoopNodeNumber], localNodeNumber);
			 pathNumArray[localNodeNumber] = pathNum; }
	       }
	  }
	  if (doMinimalWorkHere)
	       goto mustSplit2;

	  for (int i=overlapMatchIndexHold-1; i>=0; i--) {
	       std::map<unitigLocStruct, int>::iterator it;
	       int isGood = 0;
	       tempULS.ori = evenReadMatchStructs[i].ori;
	       if (evenReadMatchStructs[i].ori == 'F')
		    tempULS.frontEdgeOffset = startValue - evenReadMatchStructs[i].bhg;
	       else
		    tempULS.frontEdgeOffset = startValue + evenReadMatchStructs[i].ahg;
	       tempULS.unitig2 = evenReadMatchStructs[i].kUnitigNumber;
	       it = nodeToIndexMap.find (tempULS);
	       if (it == nodeToIndexMap.end())
		    break;
	       int nodeNum = it->second;
	       if (pathNumArray[nodeNum] == 0)
		    break;
	       isGood = 0;
	       for (int j=revStartIndices[localNodeNumberHold]; j<revStartIndices[localNodeNumberHold]+revNumIndices[localNodeNumberHold]; j++) {
		    localNodeNumber = revConnections[j].node1;
		    if (localNodeNumber == nodeNum) {
			 pathNumArray[nodeNum] = pathNum;
			 localNodeNumberHold = i;
			 isGood = 1;
			 break; }
	       }
	       if (! isGood)
		    break;
	       localNodeNumberHold = nodeNum;
	  }
	       
	  nodeIntArray.push(localNodeNumberHold);
	  while (! nodeIntArray.empty()) {
	       int localLoopNodeNumber = nodeIntArray.top();
	       nodeIntArray.pop();
//	       fprintf (stderr, "localLoopNodeNumber = %d, start index = %d, num indices = %d\n", (int) localLoopNodeNumber, (int) revStartIndices[localLoopNodeNumber], (int) revNumIndices[localLoopNodeNumber]);
	       for (int j=revStartIndices[localLoopNodeNumber]; j<revStartIndices[localLoopNodeNumber]+revNumIndices[localLoopNodeNumber]; j++) {
		    localNodeNumber = revConnections[j].node1;
//		    fprintf (stderr, "pathNum = %d, localPathVal = %d\n", (int) pathNum, (int) pathNumArray[localNodeNumber]);
		    if (pathNumArray[localNodeNumber] < pathNum) {
			 nodeIntArray.push (localNodeNumber);
//			 fprintf (stderr, "localLoopNodeNumber = %d, revStartIndices = %d, numRevIndices = %d, Pushing %d\n", localLoopNodeNumber, revStartIndices[localLoopNodeNumber], revNumIndices[localLoopNodeNumber], localNodeNumber);
			 pathNumArray[localNodeNumber] = pathNum; }
	       }
	  }
	  lastGoodNodeNumber = -1;
	  for (it1=startingNodes.begin(); it1!= startingNodes.end(); it1++) {
	       tempULS.unitig2 = it1->unitig2;
	       tempULS.frontEdgeOffset = it1->frontEdgeOffset;
	       tempULS.ori = it1->ori;
	       startingNodeNumber = nodeToIndexMap[tempULS]; }
	  if (nodeIntArray.size() > 0)
	       fprintf (stderr, "ERROR in nodeIntArray: size should be 0\n");
	  unitigNodeNumbersForPath.clear();
	  nodeIntArray.push (startingNodeNumber);
	  while (! nodeIntArray.empty()) {
	       int localLoopNodeNumber = nodeIntArray.top();
	       unitigNodeNumbersForPath.push_back(localLoopNodeNumber);
	       nodeIntArray.pop();
//             fprintf (stderr, "localLoopNodeNumber = %d, start index = %d, num indices = %d\n", (int) localLoopNodeNumber, (int) fwdStartIndices[localLoopNodeNumber], (int) fwdNumIndices[localLoopNodeNumber]);
	       for (int j=fwdStartIndices[localLoopNodeNumber]; j<fwdStartIndices[localLoopNodeNumber]+fwdNumIndices[localLoopNodeNumber]; j++) {
                    localNodeNumber = fwdConnections[j].node2;
//                  fprintf (stderr, "pathNum = %d, localPathVal = %d\n", (int) pathNum, (int) pathNumArray[localNodeNumber]);
                    if (pathNumArray[localNodeNumber] == pathNum) {
                         nodeIntArray.push (localNodeNumber);
//                       fprintf (stderr, "localLoopNodeNumber = %d, fwdStartIndices = %d, numFwdIndices = %d, Pushing %d\n", localLoopNodeNumber, fwdStartIndices[localLoopNodeNumber], fwdNumIndices[localLoopNodeNumber], localNodeNumber);
                         pathNumArray[localNodeNumber] = pathNum; }
               }
	       if (nodeIntArray.size() > 1) {
		    lastGoodNodeNumber = localLoopNodeNumber;
		    while (! nodeIntArray.empty())
			 nodeIntArray.pop();
		    break; }
          }
	  if (lastGoodNodeNumber < 0) {
	       for (numUnitigPathPrintRecsOnPath=0; numUnitigPathPrintRecsOnPath<(int) unitigNodeNumbersForPath.size(); numUnitigPathPrintRecsOnPath++) {
		    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].unitig1 = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].unitig2;
		    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].frontEdgeOffset = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].frontEdgeOffset;
		    if (numUnitigPathPrintRecsOnPath > 0)
			 augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = 1;
		    else
			 augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = 0;
		    if (numUnitigPathPrintRecsOnPath<(int)unitigNodeNumbersForPath.size() - 1)
			 augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = 1;
		    else
			 augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = 0;
//		    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].beginOffset = ???;
//		    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].endOffset = ???;
		    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].ori = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].ori;
	       }
	       generateSuperReadPlacementLinesForJoinedMates();
	       approxNumPaths = 1;
	       ++numJoinableAfterRead1Analysis;
#ifdef KILL120102
	       fprintf (stderr, "We have successfully joined the mates!\n");
	       for (unsigned int i=0; i<unitigNodeNumbersForPath.size(); i++) {
		    fprintf (stderr, "%d %c ", nodeArray[unitigNodeNumbersForPath[i]].unitig2, nodeArray[unitigNodeNumbersForPath[i]].ori);
		    fprintf (stderr, "\n"); }
#endif
	       goto afterSuperRead; }
#ifdef KILL120102
	  fprintf (stderr, "The node before the split is uni = %d offset = %d ori = %c\n", nodeArray[lastGoodNodeNumber].unitig2, nodeArray[lastGoodNodeNumber].frontEdgeOffset, nodeArray[lastGoodNodeNumber].ori);
#endif	  
	  
     mustSplit2:
	  if (oddReadMatchStructs.size() == 1)
	       goto afterSuperRead;
	  // startValue is the number of bases in the super-read beyond
	  // the last base in the read
          if (oddReadMatchStructs[0].ori == 'F')
               startValue = oddReadMatchStructs[0].ahg;
          else
               startValue = - oddReadMatchStructs[0].bhg;
//          fprintf (stderr, "Starting nodes:\n");
//          for (it1=startingNodes.begin(); it1!= startingNodes.end(); it1++)
//               fprintf (stderr, "%d %d %c\n", it1->unitig2, it1->frontEdgeOffset, it1->ori);
//          fprintf (stderr, "Ending nodes:\n");
//          for (it1=endingNodes.begin(); it1!= endingNodes.end(); it1++)
//               fprintf (stderr, "%d %d %c\n", it1->unitig2, it1->frontEdgeOffset, it1->ori);

	  // Here we are measuring the distance of the frontEdgeOffset
	  // from the end of the super-read (using a positive distance)
          for (int i=oddReadMatchStructs.size()-1; i>=0; i--) {
	       int distFromEndOfSuperRead = 0;
	       if (oddReadMatchStructs[i].ori == 'F')
		    abbULS1.ori = 'R';
	       else
		    abbULS1.ori = 'F';
	       if (oddReadMatchStructs[i].ori == 'F')
		    distFromEndOfSuperRead = startValue - oddReadMatchStructs[i].ahg;
	       else
		    distFromEndOfSuperRead = startValue + oddReadMatchStructs[i].bhg;
	       // We now have to modify the offset so as to be in relation to the
	       // left end of the super-read. We must consider all possible
	       // distances to do this.
	       numPossibleLengths = 0;
	       for (it1=endingNodes.begin(); it1!=endingNodes.end(); it1++) {
		    abbULS1.frontEdgeOffset = it1->frontEdgeOffset - distFromEndOfSuperRead;
		    elementIndex = treeFindElement (treeArr + oddReadMatchStructs[i].kUnitigNumber, (char *) &abbULS1);
		    if (elementIndex == TREE_NIL)
			 continue;
		    tempULS.unitig2 = oddReadMatchStructs[i].kUnitigNumber;
		    tempULS.frontEdgeOffset = abbULS1.frontEdgeOffset;
		    tempULS.ori = abbULS1.ori;
		    it2 = nodeToIndexMap.find (tempULS);
		    if (it2 == nodeToIndexMap.end())
			 continue;
		    if (pathNumArray[it2->second] != pathNum)
			 continue;
		    elementIndex2 = elementIndex;
		    localSuperReadLength = it1->frontEdgeOffset;
		    ++numPossibleLengths;
	       }
	       if (numPossibleLengths != 1) {
#ifdef KILL120102
		    if (numPossibleLengths == 0)
			 fprintf (stderr, "No good read2 unitig match %s %d %d %d FAIL %d %d\n", readNameSpace, treeSize, i, abbULS1.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax);
		    else
			 fprintf (stderr, "Too many good read2 unitig matches. Last is %s %d %d %d FAIL %d %d\n", readNameSpace, treeSize, i, abbULS1.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax);
#endif
		    continue;
	       }
               char *vptr;
               setTreeValPtr (vptr, treeArr+oddReadMatchStructs[i].kUnitigNumber, elementIndex2);
               // If we get here we have a unitig on the path
#ifdef KILL120102
               fprintf (stderr, "%s %d %d %d SUCCESS %d %d %d\n", readNameSpace, treeSize, i, abbULS1.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax, pathNum);
#endif
               if (abbULS1.frontEdgeOffset <= splitJoinWindowMin)
		    continue;
               // Checking for useless result
               if (abbULS1.frontEdgeOffset >= splitJoinWindowMax)
                    goto afterSuperRead;
               // If we get here we've passed all the tests and we can use it
               localUnitigNumber = oddReadMatchStructs[i].kUnitigNumber;
	       localFrontEdgeOffset = abbULS1.frontEdgeOffset;
               overlapMatchIndexHold = i;
	       break;
          }
          if (elementIndex == TREE_NIL)
               goto afterSuperRead;
	       
          tempULS.unitig2 = localUnitigNumber;
          tempULS.ori = abbULS1.ori;
          tempULS.frontEdgeOffset = localFrontEdgeOffset;
#ifdef KILL120102
          fprintf (stderr, "unitig2 = %d, frontEdgeOffset = %d, ori = %c\n", tempULS.unitig2, tempULS.frontEdgeOffset, tempULS.ori);
#endif
          localNodeNumber = localNodeNumberHold = nodeToIndexMap[tempULS];
#ifdef KILL120102
          fprintf (stderr, "localNodeNumber = %d\n", localNodeNumber);
#endif
          ++pathNum;
          pathNumArray[localNodeNumber] = pathNum;
          nodeIntArray.push(localNodeNumber);
          while (! nodeIntArray.empty()) {
               int localLoopNodeNumber = nodeIntArray.top();
               nodeIntArray.pop();
//             fprintf (stderr, "localLoopNodeNumber = %d, start index = %d, num indices = %d\n", (int) localLoopNodeNumber, (int) revStartIndices[localLoopNodeNumber], (int) revNumIndices[localLoopNodeNumber]);
               for (int j=revStartIndices[localLoopNodeNumber]; j<revStartIndices[localLoopNodeNumber]+revNumIndices[localLoopNodeNumber]; j++) {
                    localNodeNumber = revConnections[j].node1;
//                  fprintf (stderr, "pathNum = %d, localPathVal = %d\n", (int) pathNum, (int) pathNumArray[localNodeNumber]);
                    if (pathNumArray[localNodeNumber] == pathNum-1) {
                         nodeIntArray.push (localNodeNumber);
//                       fprintf (stderr, "localLoopNodeNumber = %d, revStartIndices = %d, numRevIndices = %d, Pushing %d\n", localLoopNodeNumber, revStartIndices[localLoopNodeNumber], revNumIndices[localLoopNodeNumber], localNodeNumber);
                         pathNumArray[localNodeNumber] = pathNum; }
               }
          }

          for (int i=overlapMatchIndexHold-1; i>=0; i--) {
               std::map<unitigLocStruct, int>::iterator it;
               int isGood = 0;
	       if (oddReadMatchStructs[i].ori == 'F')
		    tempULS.ori = 'R';
	       else
		    tempULS.ori = 'F';
               if (oddReadMatchStructs[i].ori == 'F')
                    distFromEndOfSuperRead = startValue - oddReadMatchStructs[i].ahg;
               else
		    distFromEndOfSuperRead = startValue - oddReadMatchStructs[i].bhg;
	       tempULS.frontEdgeOffset = localSuperReadLength - distFromEndOfSuperRead;
               tempULS.unitig2 = oddReadMatchStructs[i].kUnitigNumber;
               it = nodeToIndexMap.find (tempULS);
               if (it == nodeToIndexMap.end())
                    break;
               int nodeNum = it->second;
//               if (pathNumArray[nodeNum] == 0)
	       if (pathNumArray[nodeNum] < pathNum-1)
                    break;
               isGood = 0;
               for (int j=fwdStartIndices[localNodeNumberHold]; j<fwdStartIndices[localNodeNumberHold]+fwdNumIndices[localNodeNumberHold]; j++) {
                    localNodeNumber = fwdConnections[j].node2;
                    if (localNodeNumber == nodeNum) {
                         pathNumArray[nodeNum] = pathNum;
                         localNodeNumberHold = i;
                         isGood = 1;
                         break; }
               }
               if (! isGood)
                    break;
               localNodeNumberHold = nodeNum;
          }

          nodeIntArray.push(localNodeNumberHold);

//        fprintf (stderr, "fwdConnections\n");
//        for (unsigned int j=0; j<fwdConnections.size(); j++)
//             fprintf (stderr, "%d %d\n", fwdConnections[j].node1, fwdConnections[j].node2);
          while (! nodeIntArray.empty()) {
               int localLoopNodeNumber = nodeIntArray.top();
               nodeIntArray.pop();
//             fprintf (stderr, "localLoopNodeNumber = %d, start index = %d, num indices = %d\n", (int) localLoopNodeNumber, (int) fwdStartIndices[localLoopNodeNumber], (int) fwdNumIndices[localLoopNodeNumber]);
               for (int j=fwdStartIndices[localLoopNodeNumber]; j<fwdStartIndices[localLoopNodeNumber]+fwdNumIndices[localLoopNodeNumber]; j++) {
                    localNodeNumber = fwdConnections[j].node2;
//                  fprintf (stderr, "pathNum = %d, localPathVal = %d\n", (int) pathNum, (int) pathNumArray[localNodeNumber]);
//		    if (pathNumArray[localNodeNumber] < pathNum) {
                    if (pathNumArray[localNodeNumber] == pathNum-1) {
                         nodeIntArray.push (localNodeNumber);
//                       fprintf (stderr, "localLoopNodeNumber = %d, fwdStartIndices = %d, numFwdIndices = %d, Pushing %d\n", localLoopNodeNumber, fwdStartIndices[localLoopNodeNumber], fwdNumIndices[localLoopNodeNumber], localNodeNumber);
                         pathNumArray[localNodeNumber] = pathNum; }
               }
          }
               
	  // We now analyze the path we have to see if it's unique
          lastGoodNodeNumber = -1;
//          for (it1=startingNodes.begin(); it1!= startingNodes.end(); it1++) {
//               tempULS.unitig2 = it1->unitig2;
//               tempULS.frontEdgeOffset = it1->frontEdgeOffset;
//               tempULS.ori = it1->ori;
//               startingNodeNumber = nodeToIndexMap[tempULS]; }
          if (nodeIntArray.size() > 0)
               fprintf (stderr, "ERROR in nodeIntArray: size should be 0\n");
          unitigNodeNumbersForPath.clear();
          nodeIntArray.push (startingNodeNumber);
          while (! nodeIntArray.empty()) {
               int localLoopNodeNumber = nodeIntArray.top();
               unitigNodeNumbersForPath.push_back(localLoopNodeNumber);
               nodeIntArray.pop();
//             fprintf (stderr, "localLoopNodeNumber = %d, start index = %d, num indices = %d\n", (int) localLoopNodeNumber, (int) fwdStartIndices[localLoopNodeNumber], (int) fwdNumIndices[localLoopNodeNumber]);
               for (int j=fwdStartIndices[localLoopNodeNumber]; j<fwdStartIndices[localLoopNodeNumber]+fwdNumIndices[localLoopNodeNumber]; j++) {
                    localNodeNumber = fwdConnections[j].node2;
//                  fprintf (stderr, "pathNum = %d, localPathVal = %d\n", (int) pathNum, (int) pathNumArray[localNodeNumber]);
                    if (pathNumArray[localNodeNumber] == pathNum) {
                         nodeIntArray.push (localNodeNumber);
//                       fprintf (stderr, "localLoopNodeNumber = %d, fwdStartIndices = %d, numFwdIndices = %d, Pushing %d\n", localLoopNodeNumber, fwdStartIndices[localLoopNodeNumber], fwdNumIndices[localLoopNodeNumber], localNodeNumber);
                         pathNumArray[localNodeNumber] = pathNum; }
               }
               if (nodeIntArray.size() > 1) {
                    lastGoodNodeNumber = localLoopNodeNumber;
		    ++numJoinableUnresolvedAtEnd;
                    while (! nodeIntArray.empty())
                         nodeIntArray.pop();
                    break; }
          }
          if (lastGoodNodeNumber < 0) {
               for (numUnitigPathPrintRecsOnPath=0; numUnitigPathPrintRecsOnPath<(int) unitigNodeNumbersForPath.size(); numUnitigPathPrintRecsOnPath++) {
                    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].unitig1 = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].unitig2;
                    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].frontEdgeOffset = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].frontEdgeOffset;
                    if (numUnitigPathPrintRecsOnPath > 0)
                         augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = 1;
                    else
                         augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = 0;
                    if (numUnitigPathPrintRecsOnPath<(int)unitigNodeNumbersForPath.size() - 1)
                         augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = 1;
                    else
                         augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = 0;
//                  augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].beginOffset = ???;
//                  augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].endOffset = ???;
                    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].ori = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].ori;
               }
               generateSuperReadPlacementLinesForJoinedMates();
               approxNumPaths = 1;
	       ++numJoinableAfterBothReadAnalysis;
#ifdef KILL120102
	       fprintf (stderr, "We have successfully joined the mates!\n");
	       for (unsigned int i=0; i<unitigNodeNumbersForPath.size(); i++) {
		    fprintf (stderr, "%d %c ", nodeArray[unitigNodeNumbersForPath[i]].unitig2, nodeArray[unitigNodeNumbersForPath[i]].ori);
                    fprintf (stderr, "\n"); }
#endif
               goto afterSuperRead; }
#ifdef KILL120102
          fprintf (stderr, "The node before the split is uni = %d offset = %d ori = %c\n", nodeArray[lastGoodNodeNumber].unitig2, nodeArray[lastGoodNodeNumber].frontEdgeOffset, nodeArray[lastGoodNodeNumber].ori);
#endif
          
     afterSuperRead:
	  // Cleaning up the data structures
	  for (int j = 0; j < numTreesUsed; j++)
	       treeArr[treeReinitList[j]].root = TREE_NIL;
	  numTreesUsed = 0;
	  dataArr.arraySize = 0;
	  
#ifdef KILLED111115
	  printf ("Approx num paths returned = %d\n", approxNumPaths);
#endif
	  
	  // Doing the output (if possible)
	  
	  if (approxNumPaths == 1) {
	       fputs (outputString, outputFile);
	       return; }
#ifdef KILL120103
	  if (stderrOutputString[0] != 0)
	       fputs (stderrOutputString, stderr);
#endif
	  if (approxNumPaths < 1)
	       goto outputTheReadsIndividually;
	  // Now we move up the unitig on read1 and try again
	  

     outputTheReadsIndividually:
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold-1);
	  findSingleReadSuperReads(readNameSpace);
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
	  findSingleReadSuperReads(readNameSpace);
     }
     return;
}

bool firstNodeSort (struct nodePair val1, struct nodePair val2)
{
     return (val1.node1 < val2.node1);
}
 
bool secondNodeSort (struct nodePair val1, struct nodePair val2)
{
     return (val1.node2 < val2.node2);
}

