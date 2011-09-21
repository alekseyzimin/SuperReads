#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <set>
using namespace std;

char *line;
char *flds[3000];
char *newSuperReadNameInfoSpace;
long long newSuperReadNameInfoSpaceNeeded;
int *kUnitigLengths;

// The following returns the number of kUnitigs
int loadKUnitigs (char *filename);
int getNumReads (char *filename);
int getFldsFromLine (char *cptrHold, char **returnArrayLoc);
FILE *Fopen (const char *fn, const char *mode);

#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

int main (int argc, char **argv)
{
     int ahg, ahgHold=0, bhg, bhgHold=0;
     int firstBaseOffset, i, isFirstLine, kUniLen;
     int kUniNum, readName=0, superReadLength;
     char *kUnitigLengthsFile=argv[3], *reducedKUnitigMatchesFile=argv[1], *superReadGroupFile=argv[2];
     char *superReadsEliminatedDueToBadKUnitigMatchFile=NULL;
     char readOriInSuperRead, newSuperReadNameInfoStr[50000];
     char searchStr[100], superReadName[50000];
     int indexValue=0;
     char *newSuperReadNameInfoFlds[3000];
     unsigned char wasFound, *wasReversed;
     int kUniBegins[1000], kUniEnds[1000], superReadUnis[1000];
     FILE *infile;
     set<string> isBadSuperRead;
     int numKUnitigs, numReads;
     char *cptr, *cptr2, *cptr3;
     int numSuperReadInfoFlds=0, numKUnisInSuperRead=0;
     char **newSuperReadNameInfo;
     int numNewSuperReadNameInfoFlds=0;

     // Must allocate line, kUnitigLengths, wasReversed
     mallocOrDie (line, 1000000, char);
     
     if (argc > 4) {
	  superReadsEliminatedDueToBadKUnitigMatchFile = argv[4];
	  infile = Fopen (superReadsEliminatedDueToBadKUnitigMatchFile, "r");
	  while (fgets (line, 1000000, infile)) {
	       cptr = line;
	       if (! isdigit (*cptr))
		    continue;
	       while (*cptr) {
		    if (isspace(*cptr))
			 break;
		    ++cptr;
	       }
	       *cptr = 0;
	       isBadSuperRead.insert (string (line));
	  }
	  fclose (infile);
     }

//     printf ("Got to 2\n"); fflush (stdout);

     numKUnitigs = loadKUnitigs (kUnitigLengthsFile);

//     printf ("Got to 3\n"); fflush (stdout);

     numReads = getNumReads (superReadGroupFile);

//     printf ("Got to 4\n"); fflush (stdout);

     mallocOrDie (wasReversed, numReads, unsigned char);
     mallocOrDie (newSuperReadNameInfo, numReads, char *);
     mallocOrDie (newSuperReadNameInfoSpace, newSuperReadNameInfoSpaceNeeded, char);

//     printf ("Got to 5\n"); fflush (stdout);

     cptr3 = newSuperReadNameInfoSpace;
     infile = Fopen (superReadGroupFile, "r");
     while (fgets (line, 1000000, infile)) {
	  cptr = strstr (line, "readName = ");
	  cptr += strlen ("readName = ");
//	  printf ("Got to 11\n"); fflush (stdout);
	  readName = atoi (cptr);
//	  printf ("readName = %d\n", readName);
	  while (*cptr != ':')
	       ++cptr;
	  ++cptr;
	  while (isspace(*cptr))
	       ++cptr;
	  // Now we truncate the line after the last non-whitespace character
//	  printf ("Got to 12\n"); fflush (stdout);
	  cptr2 = cptr + strlen (cptr);
	  --cptr2;
	  while (isspace(*cptr2))
	       --cptr2;
	  cptr2[1] = 0;
	  fgets (newSuperReadNameInfoStr, 50000, infile);
	  newSuperReadNameInfoStr[strlen(newSuperReadNameInfoStr)-1] = 0;
//	  printf ("Got to 13\n"); fflush (stdout);
	  wasReversed[readName] = 0;
//	  printf ("cptr = '%s', newSuperReadNameInfoStr = '%s'\n", cptr, newSuperReadNameInfoStr);
	  if (strcmp (cptr, newSuperReadNameInfoStr) != 0)
	       wasReversed[readName] = 1;
	  newSuperReadNameInfo[readName] = cptr3;
	  strcpy (cptr3, newSuperReadNameInfoStr);
	  cptr3 += (strlen(cptr3)+1);
	  numSuperReadInfoFlds = getFldsFromLine (newSuperReadNameInfoStr, newSuperReadNameInfoFlds);
	  sprintf (superReadName, "%s%c", newSuperReadNameInfoFlds[0], newSuperReadNameInfoFlds[1][0]);
	  cptr = superReadName + strlen (superReadName);
//	  printf ("Got to 15\n"); fflush (stdout);
	  for (i=2; i<numSuperReadInfoFlds; i+=3) {
	       sprintf (cptr, "_%s_%s%c", newSuperReadNameInfoFlds[i], newSuperReadNameInfoFlds[i+1], newSuperReadNameInfoFlds[i+2][0]);
	       cptr += strlen (cptr);
	  }
//	  printf ("Got to 16\n"); fflush (stdout);
	  if (isBadSuperRead.find (string(superReadName)) != isBadSuperRead.end())
	       continue;
     }
     fclose (infile);

//     printf ("Got to 6\n"); fflush (stdout);

     infile = Fopen (reducedKUnitigMatchesFile, "r");
     isFirstLine = 1;
     wasFound = 0;
     while (fgets (line, 1000000, infile)) {
//	  printf ("Got to 7\n"); fflush (stdout);
	  if (strstr (line, "readNum") != NULL) {
	       isFirstLine = 1;
	       wasFound = 0;
	       continue; }
//	  printf ("Got to 8\n"); fflush (stdout);
	  for (cptr=line; (*cptr) && (isspace(*cptr)); cptr++) {
	  }
	  if (*cptr == 0)
	       goto processTheRead;
//	  printf ("Got to 9\n"); fflush (stdout);
	  getFldsFromLine (line, flds);
	  kUniNum = atoi (flds[6]);
	  kUniLen = atoi (flds[4]);
	  ahg = atoi (flds[0])-1;
	  bhg = kUniLen - atoi (flds[1]);
	  if (atoi(flds[2]) < atoi (flds[3])) {
	       ahg -= (atoi(flds[2])-1);
	       bhg = (atoi(flds[5]) - atoi(flds[3])) - bhg;
	  }
	  else {
	       ahg -= (atoi(flds[5])-atoi(flds[2]));
	       bhg = (atoi(flds[3])-1) - bhg;
	  }
//	  printf ("Got to 10\n"); fflush (stdout);
	  if (isFirstLine) {
	       isFirstLine = 0;
	       ahgHold = bhgHold = 1000000000; // NEWTEST
	       readName = atoi(flds[7]);
	       if (newSuperReadNameInfo[readName] == NULL)
		    newSuperReadNameInfo[readName] = (char *) "";
//	       printf ("Got to 100, readName = %d, newSuperReadNameInfo address = %p\n", readName, newSuperReadNameInfo[readName]); fflush (stdout);
	       numNewSuperReadNameInfoFlds = getFldsFromLine (newSuperReadNameInfo[readName], newSuperReadNameInfoFlds);
	       if (numNewSuperReadNameInfoFlds == 0)
		    continue;
//	       printf ("Got to 101\n"); fflush (stdout);
	       superReadLength = kUnitigLengths[atoi(newSuperReadNameInfoFlds[0])];
//	       printf ("Got to 102\n"); fflush (stdout);
	       cptr = superReadName;
//	       cptr2 = superReadUniStr;
	       sprintf (cptr, "%s%c", newSuperReadNameInfoFlds[0], newSuperReadNameInfoFlds[1][0]);
//	       printf ("Got to 103\n"); fflush (stdout);
	       cptr += strlen(cptr);
//	       sprintf (cptr2, " %s ", newSuperReadNameInfoFlds[0]);
//	       cptr2 += strlen (cptr2);
	       numKUnisInSuperRead = 0;
	       kUniBegins[numKUnisInSuperRead] = 0;
//	       printf ("Got to 104, superReadLength = %d\n", superReadLength); fflush (stdout);
	       kUniEnds[numKUnisInSuperRead] = superReadLength;
//	       if (readName == 797570) 
//		    printf ("readName = %d kUniBegins[0] = %d kUniEnds[0] = %d\n", readName, kUniBegins[0], kUniEnds[0]), fflush (stdout);
	       superReadUnis[numKUnisInSuperRead] = atoi (newSuperReadNameInfoFlds[0]);
	       ++numKUnisInSuperRead;
//	       printf ("Got to 105\n"); fflush (stdout);
	       for (i=2; i<numNewSuperReadNameInfoFlds; i+=3) {
//		    printf ("Got to 1051, i = %d, numNewSuperReadNameInfoFlds = %d, vals = '%s', '%s', and '%s'\n", i, numNewSuperReadNameInfoFlds, newSuperReadNameInfoFlds[i], newSuperReadNameInfoFlds[i+1], newSuperReadNameInfoFlds[i+2]); fflush (stdout);
		    sprintf (cptr, "_%s_%s%c", newSuperReadNameInfoFlds[i], newSuperReadNameInfoFlds[i+1], newSuperReadNameInfoFlds[i+2][0]);
		    cptr += strlen (cptr);
//	       printf ("Got to 1052\n"); fflush (stdout);
		    kUniBegins[numKUnisInSuperRead] = superReadLength - atoi (newSuperReadNameInfoFlds[i]);
//	       printf ("Got to 1053\n"); fflush (stdout);
		    superReadLength += (kUnitigLengths[atoi(newSuperReadNameInfoFlds[i+1])] - atoi(newSuperReadNameInfoFlds[i]));
//	       printf ("Got to 1054\n"); fflush (stdout);
		    superReadUnis[numKUnisInSuperRead] = atoi (newSuperReadNameInfoFlds[i+1]);
//	       printf ("Got to 1055\n"); fflush (stdout);
//		    sprintf (cptr2, "%s ", newSuperReadNameInfoFlds[i+1]);
		    kUniEnds[numKUnisInSuperRead] = superReadLength;
//	       printf ("Got to 1056\n"); fflush (stdout);
		    ++numKUnisInSuperRead;
	       }
//	       printf ("Got to 106\n"); fflush (stdout);
	  }
//	  printf ("Got to 11\n"); fflush (stdout);
	  if (numNewSuperReadNameInfoFlds == 0)
	       continue;
//	  if (readName == 797570) printf ("numKUnisInSuperRead = %d\n", numKUnisInSuperRead);
	  if (! wasFound) {
	       sprintf (searchStr, " %d ", kUniNum);
	       indexValue = -1;
	       for (i=0; i<numKUnisInSuperRead; i++)
		    if (superReadUnis[i] == kUniNum) {
			 indexValue = i;
			 break;
		    }
	       if (indexValue >= 0) {
		    ahgHold = ahg;
		    bhgHold = bhg;
		    if (wasReversed[readName]) {
			 for (i=numKUnisInSuperRead-1; i>=0; i--)
			      if (superReadUnis[i] == kUniNum) {
				   indexValue = i;
				   break; }
		    }
		    wasFound = 1;
	       }
	  }
	  continue;
	  
     processTheRead:
	  if (numNewSuperReadNameInfoFlds == 0)
	       continue;
//	  printf ("Got to 12\n"); fflush (stdout);
	  if (isBadSuperRead.find (string(superReadName)) != isBadSuperRead.end())
	       continue;
//	  printf ("Got to 13\n"); fflush (stdout);
//	  if (readName == 797570) {
//	       printf ("readName = %d, wasReversed = %d, indexValue = %d, kUniBegin = %d, kUniEnd = %d, ahgHold = %d, bhgHold = %d\n", readName, wasReversed[readName], indexValue, kUniBegins[indexValue], kUniEnds[indexValue], ahgHold, bhgHold); fflush (stdout); }
	  // The following forces a match that doesn't exist so that the
	  // program doesn't blow up
	  if (indexValue < 0) indexValue = 0;
	  if (wasReversed[readName] == 0) {
//	  printf ("Got to 141\n"); fflush (stdout);
	       firstBaseOffset = kUniBegins[indexValue]; // Usually 0
	       if (newSuperReadNameInfoFlds[indexValue*3+1][0] == 'F')
		    firstBaseOffset += ahgHold;
	       else
		    firstBaseOffset -= bhgHold; }
	  else {
//	       printf ("Got to 142\n"); fflush (stdout);
	       firstBaseOffset = kUniEnds[indexValue]; // Usu. superReadLength
	       if (newSuperReadNameInfoFlds[indexValue*3+1][0] == 'R')
		    firstBaseOffset -= ahgHold;
	       else
		    firstBaseOffset += bhgHold; }
//	  printf ("Got to 15\n"); fflush (stdout);
	  if (wasReversed[readName])
	       readOriInSuperRead = 'R';
	  else
	       readOriInSuperRead = 'F';
	  if (strlen (superReadName) > 0)
	       printf ("%d %s %d %c\n", readName, superReadName, firstBaseOffset, readOriInSuperRead);
     }
     fclose (infile);
     
     return (0);
}

int loadKUnitigs (char *filename)
{
     FILE *infile;
     int numKUnitigs, kUnitigNumber;
     int numFlds;

     infile = Fopen (filename, "r");
     fgets (line, 1000000, infile);
     numFlds = getFldsFromLine (line, flds);
     rewind (infile);
     numKUnitigs = 0;
     if (numFlds == 1) {
	  while (fgets (line, 1000000, infile))
	       ++numKUnitigs;
     }
     else {
	  while (fgets (line, 1000000, infile)) {
	       getFldsFromLine (line, flds);
	       if (atoi (flds[0]) > numKUnitigs)
		    numKUnitigs = atoi (flds[0]);
	  }
	  ++numKUnitigs; // Before just last kUnitig number; now correct
     }
     rewind (infile);
     mallocOrDie (kUnitigLengths, numKUnitigs, int);
     kUnitigNumber = 0;
     if (numFlds == 1) {
	  while (fgets (line, 1000000, infile)) {
	       kUnitigLengths[kUnitigNumber] = atoi(line);
	       ++kUnitigNumber; }
     }
     else {
	  while (fgets (line, 1000000, infile)) {
	       getFldsFromLine (line, flds);
	       kUnitigLengths[atoi(flds[0])] = atoi(flds[1]);
	  }
     }
     fclose (infile);

     return (numKUnitigs);
}

int getNumReads (char *filename)
{
     int largestReadNum = 0, readName;
     FILE *infile;
     char *cptr;

     infile = Fopen (filename, "r");

     newSuperReadNameInfoSpaceNeeded = 0;
     while (fgets (line, 1000000, infile)) {
	  cptr = strstr (line, "readName = ");
	  cptr += strlen ("readName = ");
	  readName = atoi (cptr);
	  if (readName > largestReadNum)
	       largestReadNum = readName;
	  fgets (line, 1000000, infile);
	  newSuperReadNameInfoSpaceNeeded += strlen (line);
     }
     fclose (infile);

     return (largestReadNum+1);
}

int getFldsFromLine (char *cptrHold, char **returnArrayLoc)
{
     int numFlds=0, state = 0;
     char *cptr;

     for (cptr=cptrHold; *cptr; cptr++) {
          if (isspace (*cptr)) { state = 0;*cptr = 0; }
          else {
               if (state == 1) continue;
               returnArrayLoc[numFlds] = cptr;
               ++numFlds;
               state = 1;
          }
     }
     return (numFlds);
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

