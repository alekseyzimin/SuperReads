// This changes the read numbers back to read names
// Arg 1: The file with the table indicating the 2-character read name
//   starts to the first read number for their group
// Arg 2: The file with the read numbers starting in column 1, one per line
// Output to stdout
#include <stdio.h>
#include <ctype.h>
#include "diskBasedUnitigger.h"

int main (int argc, char **argv)
{
     int numGroups=0, groupNum, numInGroup, lowIndex, highIndex, index;
     char line[1000000];
     FILE *infile;
     char *char1s, *char2s, letter1, letter2, *cptr;
     unsigned long long *startIndex, groupStartNum, lluCurrentValue, newValue;

     infile = Fopen (argv[1], "r");
     while (fgets (line, 1000000, infile))
	  ++numGroups;
     rewind (infile);
     mallocOrDie (char1s, numGroups, char);
     mallocOrDie (char2s, numGroups, char);
     mallocOrDie (startIndex, numGroups+1, unsigned long long);
     groupNum = 0;
     while (fgets (line, 1000000, infile)) {
	  sscanf (line, "%c%c %d %llu\n", &letter1, &letter2, &numInGroup, &groupStartNum);
	  char1s[groupNum] = letter1;
	  char2s[groupNum] = letter2;
	  startIndex[groupNum] = groupStartNum;
	  fprintf (stderr, "groupNum = %d\n", groupNum);
	  startIndex[groupNum+1] = groupStartNum + numInGroup;
	  ++groupNum;
     }
     fclose (infile);

     infile = Fopen (argv[2], "r");
     while (fgets (line, 1000000, infile)) {
	  sscanf (line, "%llu", &lluCurrentValue);
	  cptr = line;
	  while (! isspace (*cptr))
	       ++cptr;
	  lowIndex = 0; highIndex = numGroups+1;
	  while (highIndex-lowIndex>1) {
	       index = (lowIndex+highIndex)/2;
	       if (lluCurrentValue < startIndex[index])
		    highIndex = index;
	       else
		    lowIndex = index; }
	  newValue = lluCurrentValue - startIndex[lowIndex];
	  printf ("%c%c%llu", char1s[lowIndex], char2s[lowIndex], newValue);
	  fputs (cptr, stdout);
     }
     fclose (infile);

     return (0);
}

