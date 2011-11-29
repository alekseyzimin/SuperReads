/* This program expects 2 args:
   1) The name of the file containing the read placements in the super-reads
   2) The error file created when creating the super-read sequences
   The file containing read placements for reads in the good super-reads
      will come out on stdout.
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <set>
#include <src/charb.hpp>
using namespace std;

charb line(1000000);

FILE *Fopen (const char *fn, const char *mode);

int main (int argc, char **argv)
{
     set<string> isBadSuperRead;
     FILE *infile;
     char *errorFilename=argv[2];
     char *readPlacementFilename = argv[1];
     char *cptr, *cptr2;
     

     infile = Fopen (errorFilename, "r");
     while (fgets (line, 1000000, infile)) {
	  if (! isdigit(line[0]))
	       continue;
	  cptr = line;
	  while (! isspace(*cptr)) ++cptr;
	  *cptr = 0;
	  isBadSuperRead.insert (string(line));
     }
     fclose (infile);

     infile = Fopen (readPlacementFilename, "r");
     while (fgets (line, 1000000, infile)) {
	  cptr = line;
	  while (! isspace (*cptr)) ++cptr;
	  while (isspace (*cptr)) ++cptr;
          // Now points to super-read name
	  cptr2 = cptr;
	  while (! isspace (*cptr2)) ++cptr2;
	  *cptr2 = 0;
	  if (isBadSuperRead.find (string(cptr)) != isBadSuperRead.end())
	       continue;
	  cptr2[0] = ' ';
	  fputs (line, stdout);
     }
     fclose (infile);
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

