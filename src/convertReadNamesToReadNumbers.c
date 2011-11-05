// Takes the input fasta file and the file with the counts for each
// read name prefix and outputs a fasta file with read numbers
// Arg 1+: input fasta file
// Last arg: input prefix counts file
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

unsigned long long firstReadNumForBin[256][256];
char line[1000000];

#define outputFauxRead(num) printf (">%llu faux\nACGT\n", num);

int main (int argc, char **argv)
{
     FILE *infile;
     int readNumber;
     char letter1, letter2, *cptr;
     // I would actually set the last to -1 to be correct, but since
     // it is unsigned and we get the same result using 1, I use that instead
     unsigned long long ullnumber, globalReadNumber, globalReadNumberHold=1;
     int i;

     memset (firstReadNumForBin, 0, 256*256*sizeof(unsigned long long));
     infile = fopen (argv[argc-1], "r");
     while (fgets (line, 1000000, infile)) {
	  sscanf (line, "%c%c %*d %llu\n", &letter1, &letter2, &ullnumber);
	  firstReadNumForBin[(int)letter1][(int)letter2] = ullnumber;
     }
     fclose (infile);

     for (i=1; i<argc-1; i++) {
	  infile = fopen (argv[i], "r");
	  while (fgets (line, 1000000, infile)) {
	       if (line[0] != '>') {
		    fputs (line, stdout);
		    continue;
	       }
	       readNumber = atoi (line+3);
	       globalReadNumber = firstReadNumForBin[(int)line[1]][(int)line[2]];
	       globalReadNumber += readNumber;
	       if (globalReadNumber != globalReadNumberHold+1) {
		    if (globalReadNumberHold % 2 == 0)
			 outputFauxRead(globalReadNumberHold+1);
		    if (globalReadNumber % 2 == 1)
			 outputFauxRead(globalReadNumber-1); }
	       globalReadNumberHold = globalReadNumber;
	       fprintf (stdout, ">%llu", globalReadNumber);
	       cptr = line+4;
	       while (! isspace(*cptr))
		    ++cptr;
	       fputs (cptr, stdout);
	  }
	  fclose (infile);
     }
     
     return (0);
}

