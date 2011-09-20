// I am converting this to c since it runs very slowly under perl
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char lines[1000][1000];
int kUnitigNums[1000], minOffsets[1000], maxOffsets[1000], overlaps[1000];
char oris[1000];

char *flds[1000];

int getFldsFromLine (char *cptr);

int main(int argc, char **argv)
{
    int inputIsFromAFile=0;
    char *inputFile;
    FILE *infile, *outfile=stdout;
    char hdrLine[1000];
    int maxTotAllowableMissingOnEnds = 2;
    int numLines=0;
    int i, numFlds;
    int readLength=0, readOffset1, readOffset2, kUnitigNum, totMissingFromEnds;
    int good;
    
    if (argc >= 2) {
	 inputIsFromAFile = 1;
	 inputFile = argv[1];
	 infile = fopen (inputFile, "r"); }
    else
	 infile = stdin;
    
    while (1) {
	 if (! (fgets (lines[numLines], 1000, infile)))
	      break;
	 if (strncmp (lines[numLines], "readNum", strlen ("readNum")) == 0) {
	      strcpy (hdrLine, lines[numLines]);
	      continue; }
	 else {
	      int hasAPrintableChar = 0;
	      for (i=0; lines[numLines][i]; i++)
		   if (isprint(lines[numLines][i])) {
			hasAPrintableChar = 1;
			break; }
	      if (hasAPrintableChar) {
		   ++numLines;
		   continue; } }
	 // Use numLines to specify the number of lines in the group
	 if (numLines > 3) {
	      numLines = 0;
	      continue; }

	 for (i=0; i<numLines; i++) {
	      numFlds = getFldsFromLine(lines[i]);
	      readOffset1 = atoi (flds[2]);
	      readOffset2 = atoi (flds[3]);
	      kUnitigNum = atoi (flds[numFlds-2]);
	      if (i == 0)
		   readLength = atoi (flds[5]);
	      kUnitigNums[i] = kUnitigNum;
	      if (readOffset1 < readOffset2) {
		   oris[i] = 'F';
		   minOffsets[i] = readOffset1;
		   maxOffsets[i] = readOffset2; }
	      else {
		   oris[i] = 'R';
		   minOffsets[i] = readOffset2;
		   maxOffsets[i] = readOffset1; } }
	 totMissingFromEnds = (minOffsets[0]-1) + (readLength - maxOffsets[numLines-1]);
	 if (totMissingFromEnds > maxTotAllowableMissingOnEnds) {
	      numLines = 0; // Set num fields to 0 before; it was an error
	      numFlds = 0; // numFlds re-set to 0
	      continue; }
	 good = 1;
	 for (i=1; i<numLines; i++) {
	      overlaps[i-1] = (maxOffsets[i-1]+1) - minOffsets[i];
	      if (overlaps[i-1] < 0)
		   good = 0; }
	 if (good) {
	      hdrLine[strlen(hdrLine)-1] = 0;
	      fputs (hdrLine, outfile);
	      fputs (" : ", outfile);
	      for (i=0; i<numLines-1; i++)
		   fprintf (outfile, "%d %c %d ", kUnitigNums[i], oris[i], overlaps[i]);
	      i = numLines-1;
	      fprintf (stdout, "%d %c\n", kUnitigNums[i], oris[i]); }
	 numLines = 0;
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

