#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int ahgLocal[1000], bhgLocal[1000], begin1sLocal[1000],
     end1sLocal[1000], begin2sLocal[1000], end2sLocal[1000],
     len1sLocal[1000], len2sLocal[1000], kUnitigsLocal[1000];
char localOris[1000];
int indicesLocal[1000], keepIndexLocal[1000];

char *flds[100];
char *line;

int getFldsFromLine (char *cptr);
int localIndexCompare (int *ptr1, int *ptr2);

#define updateLocalRecords     ahgLocal[numLocalRecords] = ahg; \
bhgLocal[numLocalRecords] = bhg; \
begin1sLocal[numLocalRecords] = atoi (flds[4]);	    \
end1sLocal[numLocalRecords] = atoi (flds[5]);	    \
begin2sLocal[numLocalRecords] = atoi (flds[6]);	    \
end2sLocal[numLocalRecords] = atoi (flds[7]);	    \
len1sLocal[numLocalRecords] = atoi (flds[8]);	    \
len2sLocal[numLocalRecords] = atoi (flds[9]);		\
kUnitigsLocal[numLocalRecords] = atoi (flds[10]); \
if (begin2sLocal[numLocalRecords] < end2sLocal[numLocalRecords]) \
localOris[numLocalRecords] = 'F'; \
else \
localOris[numLocalRecords] = 'R'
    
#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

int main (int argc, char **argv)
{
    int numFlds;
    int numLocalRecords;
    int ahg, bhg, kUnitig, i;
    char readName[100], readNameHold[100];
    int maxAllowedAlignmentDiffNotToReportKUnitigReadPair=3;
    int j, firstKeepRecord, lastKeepRecord, index1, index2, firstRecordForGroup, firstIndexForGroup;

    for (i=1; i<argc; i++) {
	 if (strcmp (argv[i], "-maxoffsetdiffnottoreport") == 0) {
	      ++i;
	      maxAllowedAlignmentDiffNotToReportKUnitigReadPair = atoi (argv[i]); }
    }

    mallocOrDie (line, 1000000, char);
    fgets (line, 1000000, stdin);
    numFlds = getFldsFromLine(line);
    sscanf (flds[3], "(%d,%d", &ahg, &bhg);
    kUnitig = atoi (flds[numFlds-2]);
    strcpy (readName, flds[numFlds-1]);
    strcpy (readNameHold, readName);
    numLocalRecords = 0;
    updateLocalRecords;
    ++numLocalRecords;
    while (1) {
	int atEof = 0;
	if (!(fgets (line, 1000000, stdin))) {
	    atEof = 1;
	    goto processTheRead; }
	numFlds = getFldsFromLine(line);
	sscanf (flds[3], "(%d,%d", &ahg, &bhg);
	kUnitig = atoi (flds[2]);
	strcpy (readName, flds[numFlds-1]);
	if (strcmp (readName, readNameHold) != 0)
	    goto processTheRead;
	goto processTheRecord;
	
    processTheRead:
	// Keep the records you want and do analysis
	for (i=0; i<numLocalRecords; i++) {
	     indicesLocal[i] = i;
	     keepIndexLocal[i] =0; }
	     
	if (numLocalRecords == 1)
	     goto afterProcessingTheRead;
	qsort (indicesLocal, numLocalRecords, sizeof (int), (int (*)()) localIndexCompare);
	for (i=1; i<numLocalRecords; i++) {
	     if (kUnitigsLocal[indicesLocal[i-1]] == kUnitigsLocal[indicesLocal[i]]) {
		  if ((ahgLocal[indicesLocal[i-1]] == ahgLocal[indicesLocal[i]]) && (begin2sLocal[indicesLocal[i-1]] == begin2sLocal[indicesLocal[i]])) {
		       keepIndexLocal[i] = keepIndexLocal[i-1];
		       keepIndexLocal[i-1] = 0; }
		  else
		       keepIndexLocal[i-1] = keepIndexLocal[i] = 1;
	     }
	}

	firstRecordForGroup = 0;
	firstIndexForGroup = indicesLocal[firstRecordForGroup];
	firstKeepRecord = lastKeepRecord = -1;
	for (i=0; i<=numLocalRecords; i++) {
	     if ((i == numLocalRecords) || (kUnitigsLocal[indicesLocal[i]] != kUnitigsLocal[firstIndexForGroup])) {
		  if (firstKeepRecord >= 0) {
		       index1 = indicesLocal[firstKeepRecord];
		       index2 = indicesLocal[lastKeepRecord];
		       // The following leaves the keep index at 1 for uninteresting lines and sets to 2 for the others
//		       printf ("We got to 0, localOris = (%c,%c) ahgs = (%d,%d), startUni = (%d)\n", localOris[index1],localOris[index2],ahgLocal[index1],ahgLocal[index2], kUnitigsLocal[indicesLocal[firstRecordForGroup]]);
		       if ((localOris[index1] != localOris[index2]) || (abs (ahgLocal[index1]-ahgLocal[index2]) > maxAllowedAlignmentDiffNotToReportKUnitigReadPair)) {
//			    printf ("We got here\n");
			    for (j=firstRecordForGroup; j<i; j++)
				 if (keepIndexLocal[j])
				      keepIndexLocal[j] = 2;
		       }
		  }
		  firstRecordForGroup = i;
		  firstIndexForGroup = indicesLocal[firstRecordForGroup];
		  firstKeepRecord = -1; }
	     if (i == numLocalRecords)
		  break;
	     if (keepIndexLocal[i]) {
		  if (firstKeepRecord < 0)
		       firstKeepRecord = i;
		  lastKeepRecord = i; }
	}

	printf ("Read name = %s\n", readNameHold);
	for (i=0; i<numLocalRecords; i++) {
	     int index = indicesLocal[i];
//	     printf ("keep = %d, ", keepIndexLocal[i]);
	     if (keepIndexLocal[i] == 2)
		  printf ("myNucmerLine (ahg,bhg) = (%d,%d): %d %d %d %d %d %d %d %s\n", ahgLocal[index], bhgLocal[index], begin1sLocal[index], end1sLocal[index], begin2sLocal[index], end2sLocal[index], len1sLocal[index], len2sLocal[index], kUnitigsLocal[index], readNameHold);
//	     printf ("keep = %d, myNucmerLine (ahg,bhg) = (%d,%d): %d %d %d %d %d %d %d %s\n", keepIndexLocal[i], ahgLocal[index], bhgLocal[index], begin1sLocal[index], end1sLocal[index], begin2sLocal[index], end2sLocal[index], len1sLocal[index], len2sLocal[index], kUnitigsLocal[index], readNameHold);
	}
	printf ("\n");
	
    afterProcessingTheRead:
	// ...then
	if (atEof)
	     break;
	numLocalRecords = 0;
	strcpy (readNameHold, readName);
	
    processTheRecord:
	updateLocalRecords;
	++numLocalRecords;
    }

    return (0);
}
	    
/*
#!/usr/bin/perl
# Cat the input in. (It should be the file
# testOutput.nucmerLinesOnly.txt)

while ($line = <STDIN>) {
    chomp ($line);
    @flds = split (" ", $line);
    $ahgBhg = $flds[3];
    $kUnitig = $flds[-2];
    $read = $flds[-1];
    $index = "$kUnitig $read";
    if ($ahgBhg ne $ahgBhg{$index}) {
	++$count{$index};
	$lines{$index} .= "$line\n";
	$ahgBhg{$index} = $ahgBhg; } }

@keys = keys %count;
for (@keys) {
    $key = $_;
    next unless ($count{$key} > 1);
    ($kuni, $rd) = split (" ", $key);
    $rptUniLines{$kuni} .= $lines{$key}; }

@rptKUnis = keys %rptUniLines;
@rptKUnis = sort byNum @rptKUnis;
for (@rptKUnis) {
    $kuni = $_;
    print $rptUniLines{$kuni},"\n";
}

sub byNum
{
    return ($a <=> $b);
}

*/

int localIndexCompare (int *ptr1, int *ptr2)
{
     if (kUnitigsLocal[*ptr1] < kUnitigsLocal[*ptr2])
	  return (-1);
     if (kUnitigsLocal[*ptr1] > kUnitigsLocal[*ptr2])
	  return (1);
     if (localOris[*ptr1] < localOris[*ptr2])
	  return (-1);
     if (localOris[*ptr1] > localOris[*ptr2])
	  return (1);
     return (begin2sLocal[*ptr1] - begin2sLocal[*ptr2]);
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

