#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <charb.hpp>
#include <src2/extendSuperReadsForUniqueKmerNeighbors_cmdline.hpp>

struct kUnitigOverlapStruct {
     int kUnitig1;
     int kUnitig2;
     int ahg;
     int bhg;
     char relativeOri;
};

int main(int argc, char **argv)
{
     charb dir;
     cmdline_parse args;
     args.parse (argc, argv);
     strcpy (dir, args.dir_arg);
     charb fname;
     sprintf (fname, "%s/maxKUnitigNumber.txt", (char *)dir);
     FILE *infile = fopen (fname, "r");
     int maxUnitigNum;
     int numItems = fscanf (infile, "%d", &maxUnitigNum);
     if (numItems != 1)
	  exit (1);
     fclose (infile);

     unsigned char *numBeforeOverlaps, *numAfterOverlaps;
     int *beforeOverlap, *afterOverlap;
     char *beforeOverlapOri, *afterOverlapOri;
     numBeforeOverlaps = (unsigned char *) calloc (maxUnitigNum+1, sizeof(unsigned char));
     numAfterOverlaps = (unsigned char *) calloc (maxUnitigNum+1, sizeof(unsigned char));
     beforeOverlap = (int *) calloc (maxUnitigNum+1, sizeof(int));
     afterOverlap = (int *) calloc (maxUnitigNum+1, sizeof(int));
     beforeOverlapOri = (char *) calloc (maxUnitigNum+1, sizeof (char));
     afterOverlapOri = (char *) calloc (maxUnitigNum+1, sizeof (char));
     for (int i=0; i<=maxUnitigNum; ++i)
	  beforeOverlap[i] = afterOverlap[i] = -1;

     kUnitigOverlapStruct kUnitigOverlapRecord;
     sprintf (fname, "%s/overlap.overlaps", (char *)dir);
     infile = fopen (fname, "rb");
     while (fread (&kUnitigOverlapRecord, sizeof (struct kUnitigOverlapStruct), 1, infile)) {
	  int localKUnitig = kUnitigOverlapRecord.kUnitig1;
	  char localOri;
	  if (kUnitigOverlapRecord.relativeOri == 'N')
	       localOri = 'F';
	  else
	       localOri = 'R';
          if (kUnitigOverlapRecord.ahg > 0) { // Forward direction of read
	       if ((kUnitigOverlapRecord.kUnitig2 != afterOverlap[localKUnitig]) || (localOri != afterOverlapOri[localKUnitig])) {
		    if (numAfterOverlaps[localKUnitig] == 0) {
			 afterOverlap[localKUnitig] = kUnitigOverlapRecord.kUnitig2;
			 afterOverlapOri[localKUnitig] = localOri; }
		    if (numAfterOverlaps[localKUnitig] < 2)
			 ++numAfterOverlaps[localKUnitig]; } }
	  else { // Reverse direction of read
	       if ((kUnitigOverlapRecord.kUnitig2 != beforeOverlap[localKUnitig]) || (localOri != beforeOverlapOri[localKUnitig])) {
		    if (numBeforeOverlaps[localKUnitig] == 0) {
			 beforeOverlap[localKUnitig] = kUnitigOverlapRecord.kUnitig2;
			 beforeOverlapOri[localKUnitig] = localOri; }
		    if (numBeforeOverlaps[localKUnitig] < 2)
			 ++numBeforeOverlaps[localKUnitig]; } }
     }
     fclose (infile);

     sprintf (fname, "%s/overlap.selfOverlaps.txt", (char *)dir);
     infile = fopen (fname, "r");
     charb line(100);
     while (fgets (line, 100, infile)) {
	  int localKUnitig, localKUnitig2;
	  char localOri;
	  sscanf (line, "%d %d %c", &localKUnitig, &localKUnitig2, &localOri);
	  if (localOri == 'N')
	       localOri = 'F';
	  else
	       localOri = 'R';
	  numBeforeOverlaps[localKUnitig] = numAfterOverlaps[localKUnitig] = 2;
	  beforeOverlap[localKUnitig] = afterOverlap[localKUnitig] = -1;
	  beforeOverlapOri[localKUnitig] = afterOverlapOri[localKUnitig] = localOri;
     }
     fclose (infile);

     for (int localKUnitig=0; localKUnitig<=maxUnitigNum; ++localKUnitig) {
	  if (numBeforeOverlaps[localKUnitig] > 0) {
	       printf ("BEF %d %d%c", localKUnitig, beforeOverlap[localKUnitig], beforeOverlapOri[localKUnitig]);
	       if (numBeforeOverlaps[localKUnitig] > 1)
		    printf (" BAD");
	       printf ("\n"); }
	  if (numAfterOverlaps[localKUnitig] > 0) {
	       printf ("AFT %d %d%c", localKUnitig, afterOverlap[localKUnitig], afterOverlapOri[localKUnitig]);
	       if (numAfterOverlaps[localKUnitig] > 1)
		    printf (" BAD");
	       printf ("\n"); }
     }

     return (0);
}

