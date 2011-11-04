#include <stdio.h>
#include <stdlib.h>

struct goodPairInfoStruct {
     unsigned int pairNumber;
     int read1;
     int read2;
     short offset;
     unsigned short overlap;
};

/* We discuss the field extraElimCode in the following structure here.
   0x01: kill to begin of lower-numbered read of pairNumber1
   0x02: kill to end of lower-numbered read of pairNumber1
   0x04: kill to begin of higher-numbered read of pairNumber1
   0x08: kill to end of higher-numbered read of pairNumber1
   0x10: kill to begin of lower-numbered read of pairNumber2
   0x20: kill to end of lower-numbered read of pairNumber2
   0x40: kill to begin of higher-numbered read of pairNumber2
   0x80: kill to end of higher-numbered read of pairNumber2
 */
struct pairsToCheckStruct
{
     unsigned int pairNumber1;
     unsigned int pairNumber2;
     int read1;
     int read2;
     short offset;
     unsigned char pairGroup;
     /* The following keeps track of if we need to kill to begin or 
	end of any read. See above */
     unsigned char extraElimCode;
};

struct readIntervalStruct {
     int readNo;
     unsigned short begin;
     unsigned short end;
};

struct failingPairStruct {
     unsigned int pairNumber;
     unsigned char extraKillCode;
};

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

FILE *Popen (const char *fn, const char *mode)
{
     FILE *result;
     result = popen (fn, mode);
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

int Fwrite (const void *ptr, int sizeOfElement, int numElements, FILE *fptr)
{
     int retcode;
     retcode = fwrite (ptr, sizeOfElement, numElements, fptr);
     if (retcode != numElements)
     {
          fprintf (stderr, "Write error! Presume disk is full. Bye!\n");
          exit (2);
     }

     return (retcode);
}

int getInt (char *fname)
{
     FILE *infile;
     int tval;

     infile = Fopen (fname, "r");
     int fields_read = fscanf (infile, "%d\n", &tval);
     if(fields_read != 1) {
       fprintf(stderr, "Failed to read one int. Bye!\n");
       exit(2);
     }
     fclose (infile);
     return (tval);
}

#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

#define setTwoDimensionalMatrix(name, tempArrayName, numFirstIndex, numSecondIndex, type) mallocOrDie(tempArrayName, (numFirstIndex * numSecondIndex), type); \
mallocOrDie(name, numFirstIndex, type *); \
{ int ijkl; for (ijkl=0; ijkl<numFirstIndex; ijkl++) { name[ijkl] = &(tempArrayName[ijkl*(numSecondIndex)]); } }

