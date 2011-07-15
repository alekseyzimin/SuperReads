/* Replaces combineConsensusResultsForReadGroup.perl
   to now work with my version of read matches (was show-coords
   output before.)
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN_MATCH_LENGTH 32

int beginUnitigOffsets[1000], endUnitigOffsets[1000], beginReadOffsets[1000], endReadOffsets[1000], unitigLengths[1000], readLengths[1000], unitigNumbers[1000];
int minForLine[1000], maxForLine[1000];
int recordIndices[1000];

int indexCompare (int *iptr1, int *iptr2);

int main (int argc, char **argv)
{
     FILE *infile, *outfile;
     int readNum, recordNumForRead;
     char line[1000], *cptr;
     char readName[1000], newReadName[1000];
     int i;
     int atEof=0;

     if (argc >= 2)
	  infile = fopen (argv[1], "r");
     else
	  infile = stdin;
     if (argc >= 3)
	  outfile = fopen (argv[2], "w");
     else
	  outfile = stdout;

     readNum = 0;
     recordNumForRead = 0;
     fgets (line, 1000, infile);
     cptr = strchr (line, ':');
     ++cptr;

     sscanf (cptr, "%d %d %d %d %d %d %d %s", (beginUnitigOffsets+recordNumForRead), (endUnitigOffsets+recordNumForRead), (beginReadOffsets+recordNumForRead), (endReadOffsets+recordNumForRead), (unitigLengths+recordNumForRead), (readLengths+recordNumForRead), (unitigNumbers+recordNumForRead), newReadName);
     if (endUnitigOffsets[recordNumForRead]-beginUnitigOffsets[recordNumForRead]+1 >= MIN_MATCH_LENGTH) {
	  if (beginReadOffsets[recordNumForRead] < endReadOffsets[recordNumForRead]) {
	       minForLine[recordNumForRead] = beginReadOffsets[recordNumForRead];
	       maxForLine[recordNumForRead] = endReadOffsets[recordNumForRead]; }
	  else {
	       minForLine[recordNumForRead] = endReadOffsets[recordNumForRead];
	       maxForLine[recordNumForRead] = beginReadOffsets[recordNumForRead]; }
	  ++recordNumForRead;
     }
     
     strcpy (readName, newReadName);

     while (1) {
	  fflush (stdout);
	  if (! fgets (line, 1000, infile)) {
	       atEof = 1;
	       goto processTheRead; }
	  cptr = strchr (line, ':');
	  ++cptr;
	  sscanf (cptr, "%d %d %d %d %d %d %d %s", (beginUnitigOffsets+recordNumForRead), (endUnitigOffsets+recordNumForRead), (beginReadOffsets+recordNumForRead), (endReadOffsets+recordNumForRead), (unitigLengths+recordNumForRead), (readLengths+recordNumForRead), (unitigNumbers+recordNumForRead), newReadName);
	  if (endUnitigOffsets[recordNumForRead]-beginUnitigOffsets[recordNumForRead]+1 < MIN_MATCH_LENGTH)
	       continue;
	  if (strcmp(readName, newReadName) == 0)
	       goto addLineToTheData;
     processTheRead:
	  for (i=0; i<recordNumForRead; i++)
	       recordIndices[i] = i;
	  qsort (recordIndices, recordNumForRead, sizeof (int), (int (*)()) indexCompare);
	  fprintf (outfile, "readNum = %d; readName = %s\n", readNum, readName);
	  for (i=0; i<recordNumForRead; i++) {
	       int newIndex = recordIndices[i];
	       fprintf (outfile, "%6d %6d %6d %6d %6d %6d %9d %s\n", beginUnitigOffsets[newIndex], endUnitigOffsets[newIndex], beginReadOffsets[newIndex], endReadOffsets[newIndex], unitigLengths[newIndex], readLengths[newIndex], unitigNumbers[newIndex], readName); }
	  fprintf (outfile, "\n");
	  if (atEof)
	       break;
	  recordNumForRead = 0;
	  sscanf (cptr, "%d %d %d %d %d %d %d %s", (beginUnitigOffsets+recordNumForRead), (endUnitigOffsets+recordNumForRead), (beginReadOffsets+recordNumForRead), (endReadOffsets+recordNumForRead), (unitigLengths+recordNumForRead), (readLengths+recordNumForRead), (unitigNumbers+recordNumForRead), newReadName);
	  strcpy (readName, newReadName);
	  ++readNum;
	  if (readNum % 100000 == 0)
	       fprintf (stderr, "\r%d reads done...", readNum);
     addLineToTheData:
	  if (beginReadOffsets[recordNumForRead] < endReadOffsets[recordNumForRead]) {
	       minForLine[recordNumForRead] = beginReadOffsets[recordNumForRead];
	       maxForLine[recordNumForRead] = endReadOffsets[recordNumForRead]; }
	  else {
	       minForLine[recordNumForRead] = endReadOffsets[recordNumForRead];
	       maxForLine[recordNumForRead] = beginReadOffsets[recordNumForRead]; }
	  ++recordNumForRead; }
     fprintf (stderr, "\nFinished.\n");

     if (argc >= 2)
	  fclose (infile);
     if (argc >= 3)
	  fclose (outfile);

     return (0);
}

int indexCompare (int *iptr1, int *iptr2)
{
     if (minForLine[*iptr1] != minForLine[*iptr2])
	  return (minForLine[*iptr1] - minForLine[*iptr2]);
     return (maxForLine[*iptr2] - maxForLine[*iptr1]);
}

