#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char *flds[100];
char *line;
int *kUniLengths;
char *superReadName;
int numSuperReadPieces;

struct readPlacementDataStruct {
     int offset;
     char ori;
} readPlacementInfo[2];

int getFldsFromLine (char *cptrHold);
FILE *Fopen (const char *fn, const char *mode);
int getSuperReadLength (char *superReadName);

#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }


int main (int argc, char **argv)
{
     char *kUnitigLengthsFile = argv[1], *readPlacementFile = argv[2];
     FILE *infile;
     int maxKUniNum, kUniNum;
     long long readNum, readNumHold=-2, insertEndForRead;
     int unjoinedMateDist = -1000000000, mateDist, beginOffset, lengthFromEnd;
     int superReadLength, isFirstLine;
     int numFlds;
     
     mallocOrDie (line, 1000000, char);
     infile = Fopen (kUnitigLengthsFile, "r");
     maxKUniNum = 0; isFirstLine = 1;
     while (fgets (line, 1000000, infile)) {
	  numFlds = getFldsFromLine (line);
	  if (numFlds == 1) {
	       if (isFirstLine)
		    isFirstLine = 0;
	       else
		    ++maxKUniNum;
	       continue;
	  }
	  kUniNum = atoi (flds[0]);
	  if (kUniNum > maxKUniNum)
	       maxKUniNum = kUniNum;
     }
     rewind (infile);

     mallocOrDie (kUniLengths, maxKUniNum+1, int);

     kUniNum = 0; isFirstLine = 1;
     while (fgets (line, 1000000, infile)) {
	  numFlds = getFldsFromLine (line);
	  if (numFlds == 1) {
	       if (isFirstLine)
		    isFirstLine = 0;
	       else
		    ++kUniNum;
	       kUniLengths[kUniNum] = atoi (flds[0]);
	       continue; }
	  kUniLengths[atoi(flds[0])] = atoi (flds[1]);
     }
     fclose (infile);

     infile = Fopen (readPlacementFile, "r");
     mallocOrDie (superReadName, 50000, char);
     while (fgets (line, 1000000, infile)) {
	  getFldsFromLine (line);
	  readNum = atoll (flds[0]);
	  insertEndForRead = readNum % 2;
	  readPlacementInfo[insertEndForRead].ori = flds[3][0];
	  readPlacementInfo[insertEndForRead].offset = atoi (flds[2]);
	  if (insertEndForRead == 0) {
	       readNumHold = readNum;
	       strcpy (superReadName, flds[1]);
	       continue;
	  }
	  // Only get here if we are the second read of an insert
	  if (readNum != readNumHold+1)
	       continue;
	  // If we get here, both reads of the mate pair were placed
	  if (strcmp (superReadName, flds[1]) == 0) {
	       if ((readPlacementInfo[0].ori != readPlacementInfo[1].ori) &&
		   (readPlacementInfo[0].offset >= 0) &&
		   (readPlacementInfo[1].offset >= 0)) {
		    if (readPlacementInfo[0].ori == 'F') {
			 beginOffset = readPlacementInfo[0].offset;
			 mateDist = readPlacementInfo[1].offset - beginOffset; }
		    else {
			 beginOffset = readPlacementInfo[1].offset;
			 mateDist = readPlacementInfo[0].offset - beginOffset; }
		    // The following also sets numSuperReadPieces
		    superReadLength = getSuperReadLength (superReadName);
		    lengthFromEnd = superReadLength - beginOffset;
		    printf ("%lld %d %d %d\n", readNumHold, mateDist, lengthFromEnd, numSuperReadPieces);
	       }
	  }
	  else {
	       if ((readPlacementInfo[0].ori == 'F') && (readPlacementInfo[0].offset >= 0)) {
		    beginOffset = readPlacementInfo[0].offset;
		    superReadLength = getSuperReadLength (superReadName);
		    lengthFromEnd = superReadLength - beginOffset;
		    mateDist = unjoinedMateDist;
		    printf ("%lld %d %d %d\n", readNumHold, mateDist, lengthFromEnd, numSuperReadPieces);
	       }
	       if ((readPlacementInfo[1].ori == 'F') && (readPlacementInfo[1].offset >= 0)) {
		    beginOffset = readPlacementInfo[1].offset;
		    superReadLength = getSuperReadLength (flds[1]);
		    lengthFromEnd = superReadLength - beginOffset;
		    mateDist = unjoinedMateDist;
		    printf ("%lld %d %d %d\n", readNum, mateDist, lengthFromEnd, numSuperReadPieces);
	       }
	  }
     }

     return (0);
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

int getSuperReadLength (char *localSuperReadName)
{
     char *cptr;
     int kUnitig, superReadLength;

     cptr = localSuperReadName;
     numSuperReadPieces = 0;
     kUnitig = atoi (cptr);
     ++numSuperReadPieces;
     superReadLength = kUniLengths[kUnitig];
     ++cptr;
     while (1) {
	  while (isalnum(*cptr)) 
	       ++cptr;
	  if (*cptr == 0)
	       break;
	  // We should be at an underscore here ('_')
	  ++cptr;
	  superReadLength -= atoi(cptr);
	  while (isalnum(*cptr))
	       ++cptr;
	  // At another underscore here
	  ++cptr;
	  // ... and now at the new k-unitig
	  ++numSuperReadPieces;
	  kUnitig = atoi (cptr);
	  superReadLength += kUniLengths[kUnitig];
	  ++cptr;
     }
	  
     return (superReadLength);
}

