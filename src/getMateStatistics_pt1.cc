#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char *flds[100];
char *line;
int *kUniLengths;
char *superReadNameSpace;
int numSuperReadPieces;

struct readPlacementDataStruct {
     int offset;
     char *superReadName;
     char ori;
} *readPlacementInfo;

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
     long long maxReadNum, readNum, otherReadNum;
     long long superReadNameSpaceNeeded = 0;
     struct readPlacementDataStruct *pRPDS1, *pRPDS2;
     char *cptr;
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
     maxReadNum = 0LL;
     while (fgets (line, 1000000, infile)) {
	  getFldsFromLine (line);
	  readNum = atoll (flds[0]);
	  if (readNum > maxReadNum)
	       maxReadNum = readNum;
	  superReadNameSpaceNeeded += (strlen(flds[1]) + 1);
     }
     rewind (infile);
     mallocOrDie (superReadNameSpace, superReadNameSpaceNeeded, char);
     mallocOrDie (readPlacementInfo, maxReadNum, struct readPlacementDataStruct);
     cptr = superReadNameSpace;
     while (fgets (line, 1000000, infile)) {
	  getFldsFromLine (line);
	  readNum = atoll (flds[0]);
	  readPlacementInfo[readNum].ori = flds[3][0];
	  readPlacementInfo[readNum].offset = atoi (flds[2]);
	  strcpy (cptr, flds[1]);
	  readPlacementInfo[readNum].superReadName = cptr;
	  cptr += (strlen(cptr)+1);
     }
     fclose (infile);

     for (readNum=0; readNum<maxReadNum; readNum+=2) {
	  otherReadNum = readNum+1;
	  if (! readPlacementInfo[readNum].ori)
	       continue;
	  if (! readPlacementInfo[otherReadNum].ori)
	       continue;
	  pRPDS1 = readPlacementInfo+readNum;
	  pRPDS2 = readPlacementInfo+otherReadNum;
	  if (strcmp (pRPDS1->superReadName, pRPDS2->superReadName) == 0) {
	       if ((pRPDS1->ori != pRPDS2->ori) &&
		   (pRPDS1->offset >= 0) &&
		   (pRPDS2->offset >= 0)) {
		    if (pRPDS1->ori == 'F') {
			 beginOffset = pRPDS1->offset;
			 mateDist = pRPDS2->offset - beginOffset; }
		    else {
			 beginOffset = pRPDS2->offset;
			 mateDist = pRPDS1->offset - beginOffset; }
		    // The following also sets numSuperReadPieces
		    superReadLength = getSuperReadLength (pRPDS1->superReadName);
		    lengthFromEnd = superReadLength - beginOffset;
		    printf ("%lld %d %d %d\n", readNum, mateDist, lengthFromEnd, numSuperReadPieces);
	       }
	  }
	  else {
	       if ((pRPDS1->ori == 'F') && (pRPDS1->offset >= 0)) {
		    beginOffset = pRPDS1->offset;
		    superReadLength = getSuperReadLength (pRPDS1->superReadName);
		    lengthFromEnd = superReadLength - beginOffset;
		    mateDist = unjoinedMateDist;
		    printf ("%lld %d %d %d\n", readNum, mateDist, lengthFromEnd, numSuperReadPieces);
	       }
	       if ((pRPDS2->ori == 'F') && (pRPDS2->offset >= 0)) {
		    beginOffset = pRPDS2->offset;
		    superReadLength = getSuperReadLength (pRPDS2->superReadName);
		    lengthFromEnd = superReadLength - beginOffset;
		    mateDist = unjoinedMateDist;
		    printf ("%lld %d %d %d\n", otherReadNum, mateDist, lengthFromEnd, numSuperReadPieces);
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

int getSuperReadLength (char *superReadName)
{
     char *cptr;
     int kUnitig, superReadLength;

     cptr = superReadName;
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

