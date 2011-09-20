#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char *flds[1000];
int fldNums[100], fldNums2[100];
char *line;

int getFldsFromLine (char *cptr);
    
#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

int main (int argc, char **argv)
{
    int numFlds;
    char *cptr;
    int i, lowerValuedRecord;

    mallocOrDie (line, 1000000, char);

    while (fgets (line, 1000000, stdin)) {
	 fputs (line, stdout);
	 cptr = strrchr (line, ':');
	 ++cptr;
	 numFlds = getFldsFromLine(cptr);
	 for (i=0; i<numFlds; i+=3) {
	      fldNums[i] = atoi (flds[i]);
	      fldNums2[numFlds-2-i] = fldNums[i]; }
	 for (i=1; i<numFlds; i+=3) {
	      fldNums[i] = flds[i][0];
	      if (fldNums[i] == 'F')
		   fldNums2[numFlds-i] = 'R';
	      else
		   fldNums2[numFlds-i] = 'F'; }
	 for (i=2; i<numFlds; i+=3) {
	      fldNums[i] = atoi (flds[i]);
	      fldNums2[numFlds-i-1] = fldNums[i]; }
	 lowerValuedRecord = 0;
	 for (i=0; i<numFlds; i++) {
	      if (fldNums[i] < fldNums2[i]) {
		   lowerValuedRecord = 1;
		   break; }
	      if (fldNums[i] > fldNums2[i]) {
		   lowerValuedRecord = 2;
		   break; } }
	 if (lowerValuedRecord == 1) {
	      fprintf (stdout, "%d %c", fldNums[0], fldNums[1]);
	      for (i=2; i<numFlds; i+=3)
		   fprintf (stdout, " %d %d %c", fldNums[i], fldNums[i+1], fldNums[i+2]); }
	 else {
	      fprintf (stdout, "%d %c", fldNums2[0], fldNums2[1]);
	      for (i=2; i<numFlds; i+=3)
		   fprintf (stdout, " %d %d %c", fldNums2[i], fldNums2[i+1], fldNums2[i+2]); }
	 fprintf (stdout, "\n");
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

