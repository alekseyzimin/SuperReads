#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <vector>
#include <misc.hpp>
#include <charb.hpp>

int main (int argc, char **argv)
{
     char *fn = argv[1];
     int numNsToBreakScaffold;
     if (argc <= 2)
	  numNsToBreakScaffold = 10;
     else
	  numNsToBreakScaffold = atoi (argv[2]);
     if (numNsToBreakScaffold <= 0)
	  numNsToBreakScaffold = 1;
     int contigNum = 0;
     FILE *infile = fopen (fn, "r");
     FILE *outfile = fopen ("splitScaffoldPlacementFile.txt", "w");
     bool isFirstLine = true;
     charb line(100);
     std::vector<char *> flds;
     int numBadInARow = -1;
     charb badSequence(50);
     charb scaff(30);
     charb contigName(30);
     int curScaffOffset = 0;
     int beginContigOffsetInScaff = 0, endContigOffsetInScaff = 0;
     bool contigIsStarted = true;
     while (fgets (line, 100, infile)) {
	  if (line[0] == '>') {
	       if (! isFirstLine) {
		    fputc ('\n', stdout);
		    ++endContigOffsetInScaff;
		    fprintf (outfile, "%s %s %d %d f\n", (char *)contigName, (char *)scaff, beginContigOffsetInScaff, endContigOffsetInScaff); // Must go to file
	       }
	       else
		    isFirstLine = false;
	       contigIsStarted = false;
	       getFldsFromLine (line, flds);
	       char *cptr = flds[0]+1;
	       if (*cptr)
		    strcpy (scaff, cptr);
	       else
		    strcpy (scaff, "scaff");
	       sprintf (contigName, "%s_%d", (char *)scaff, contigNum);
	       printf (">%s\n", (char *)contigName);
	       ++contigNum;
	       numBadInARow = -1; // Negative until the first good base
	       curScaffOffset = -1;
	       continue; }
	  for (char *cptr = line; *cptr; ++cptr) {
	       if (! isspace (*cptr))
		    ++curScaffOffset;
	       switch (*cptr) {
	       case 'A': case 'a':
	       case 'C': case 'c':
	       case 'G': case 'g':
	       case 'T': case 't':
		    if (! contigIsStarted) {
			 beginContigOffsetInScaff = curScaffOffset;
			 contigIsStarted = true; }
		    if (numBadInARow > 0) {
			 if (numBadInARow >= numNsToBreakScaffold) {
			      fputc ('\n', stdout);
			      ++endContigOffsetInScaff;
			      fprintf (outfile, "%s %s %d %d f\n", (char *)contigName, (char *)scaff, beginContigOffsetInScaff, endContigOffsetInScaff); // Must go to file
			      sprintf (contigName, "%s_%d", (char *)scaff, contigNum);
			      printf (">%s\n", (char *)contigName);
			      beginContigOffsetInScaff = curScaffOffset;
			      ++contigNum; }
			 else 
			      fputs ((char *)badSequence, stdout);
			 badSequence.clear(); }
		    endContigOffsetInScaff = curScaffOffset;
		    numBadInARow = 0;
		    fputc (*cptr, stdout);
		    break;
	       default:
		    if (isspace (*cptr))
			 break;
		    if (numBadInARow < 0)
			 break;
		    ++numBadInARow;
		    if (numBadInARow < numNsToBreakScaffold)
			 strcat (badSequence, "N"); // Replacing any form of bad sequence with N's
		    break;
	       }
	  } // Matches for (char *cptr = line...
     } // Matches fgets line

     if (! isFirstLine) {
	  fputc ('\n', stdout);
	  ++endContigOffsetInScaff;
	  fprintf (outfile, "%s %s %d %d f\n", (char *)contigName, (char *)scaff, beginContigOffsetInScaff, endContigOffsetInScaff); // Must go to file
     }

     fclose (infile);
     fclose (outfile);

     return (0);
}

