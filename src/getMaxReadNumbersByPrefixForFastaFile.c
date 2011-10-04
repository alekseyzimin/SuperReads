#include <stdio.h>
#include <stdlib.h>
#include <string.h>

unsigned int numReadNumsNeededForBin[256][256];
char line[1000000];
int main (int argc, char **argv)
{
     FILE *infile;
     int readNumber, i, j;
     unsigned long long startOffset=0;

     memset (numReadNumsNeededForBin, 0, 256*256*sizeof(unsigned int));
     for (i=1; i<argc; i++) {
	  infile = fopen (argv[i], "r");
	  while (fgets (line, 1000000, infile)) {
	       if (line[0] != '>')
		    continue;
	       readNumber = atoi (line+3);
	       if (readNumber+1 > numReadNumsNeededForBin[(int)line[1]][(int)line[2]])
		    numReadNumsNeededForBin[(int)line[1]][(int)line[2]] = readNumber+1;
	  }
	  fclose (infile);
     }

     for (i=0; i<256; i++)
	  for (j=0; j<256; j++) {
	       if (numReadNumsNeededForBin[i][j] == 0)
		    continue;
	       if (numReadNumsNeededForBin[i][j] % 2 == 1)
		    ++ numReadNumsNeededForBin[i][j];
	       printf ("%c%c %d %llu\n", i, j, numReadNumsNeededForBin[i][j], startOffset);
	       startOffset += numReadNumsNeededForBin[i][j];
	  }

     return (0);
}

