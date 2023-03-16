#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>
#include <map>
#include <misc.hpp>
#include <charb.hpp>

int main (int argc, char **argv)
{
     charb line(100), querySeqNameHold, refSeqNameHold;
     charb indelOutputSequence(100), indelOutputSequenceLine(100);
     char *referenceFile = NULL, *queryFile = NULL;
     std::vector<char *>flds;
     std::map<std::string, int> len;
     int isFirstRecord = 1;
     const int FIRST_IN_INSERT = 1, SECOND_IN_INSERT = 2;
     // I haven't yet implemented indel lines, so this can wait
     bool outputOnlyPairedReads = false, noIndelLines = false;
     if (argc > 1) {
	  for (int i=1; i<argc; i++) {
	       if (strcmp (argv[i], "-mated-pairs-only") == 0)
		    outputOnlyPairedReads = true;
	       else if (strcmp (argv[i], "-no-indel-lines") == 0)
		    noIndelLines = true;
	       else if (strcmp (argv[i], "-reference-file") == 0) {
		    ++i;
		    referenceFile = argv[i]; }
	       else if (strcmp (argv[i], "-query-file") == 0) {
		    ++i;
		    queryFile = argv[i]; }
	  }
     }
     
     if (referenceFile == NULL)
	  printf ("file1 ");
     else
	  printf ("%s ", referenceFile);
     if (queryFile == NULL)
	  printf ("file2\n");
     else
	  printf ("%s\n", queryFile);
     printf ("NUCMER\n");
//     printf ("file1 file2\nNUCMER\n");
     while (fgets (line, 100, stdin)) {
	  if (line[0] == '@') {
	       getFldsFromLine (line, flds, "\t");
	       if (strcmp (flds[0], "@SQ") != 0)
		    continue;
	       if (strncmp (flds[2], "LN:", strlen("LN:")) != 0)
		    continue;
	       len[std::string(flds[1]+3)] = atoi (flds[2]+3);
	       continue; }
	  getFldsFromLine (line, flds, "\t");
	  int flagField = atoi (flds[1]);
	  int tempField = flagField % 8;
	  if (tempField >= 4) // It didn't get mapped
	       continue;
	  tempField = flagField % 512;
	  if (flagField != tempField) // Not passing quality controls
	       continue;
	  if (outputOnlyPairedReads) {
	       tempField = flagField % 4;
	       if (tempField < 2)
		    continue; }
	  int leftPos = atoi (flds[3]);
	  if (leftPos == 0)
	       continue;
	  int isRevComp = flagField % 32;
	  if (isRevComp < 16)
	       isRevComp = 0;
	  int whichEnd;
	  tempField = flagField % 256;
	  if (tempField >= 128)
	       whichEnd = SECOND_IN_INSERT;
	  else if (tempField % 128 >= 64)
	       whichEnd = FIRST_IN_INSERT;
	  int curNum = 0, currentOffset = 1, beginOffsetInQuery = 1, endOffsetInQuery = 1;
	  int currentReferenceOffset = leftPos, endOffsetInReference = leftPos;
	  bool matchHasStarted = false;
	  int lastIndelOffsetCount = 1;
	  for (char *cptr = flds[5]; *cptr; ++cptr) { // Parsing the CIGAR string
	       if (isdigit (*cptr)) {
		    curNum *= 10;
		    curNum += (*cptr - '0');
		    continue; }
	       switch (*cptr) {
	       case 'M':
	       case '=':
		    if (! matchHasStarted) {
			 matchHasStarted = true;
			 lastIndelOffsetCount = currentOffset;
			 strcpy (indelOutputSequence, "");
			 beginOffsetInQuery = currentOffset; }
	       case 'X':
		    currentOffset += curNum;
		    currentReferenceOffset += curNum;
		    if (*cptr == 'X')
			 break;
		    endOffsetInQuery = currentOffset;
		    endOffsetInReference = currentReferenceOffset;
		    break;
	       case 'I':
		    sprintf (indelOutputSequenceLine, "%d\n", -(currentOffset - lastIndelOffsetCount + 1));
		    strcat (indelOutputSequence, indelOutputSequenceLine);
		    for (int i=1; i<curNum; i++)
			 strcat (indelOutputSequence, "-1\n");
	       case 'S':
	       case 'H':
		    currentOffset += curNum;
		    if (*cptr == 'I')
			 lastIndelOffsetCount = currentOffset;
		    break;
	       case 'D':
	       case 'N':
		    sprintf (indelOutputSequenceLine, "%d\n", + (currentOffset - lastIndelOffsetCount + 1));
		    strcat (indelOutputSequence, indelOutputSequenceLine);
		    for (int i=1; i<curNum; i++)
			 strcat (indelOutputSequence, "1\n");
		    currentReferenceOffset += curNum;
		    lastIndelOffsetCount = currentOffset;
		    break;
	       }
	       curNum = 0;
	  }
	  int qBegin, qEnd, queryLength = currentOffset-1;
	  --endOffsetInReference;
	  --endOffsetInQuery;
	  if (isRevComp) {
//	       qBegin = endOffsetInQuery;
//	       qEnd = beginOffsetInQuery; }
	       qBegin = queryLength - (beginOffsetInQuery-1);
	       qEnd = (queryLength - endOffsetInQuery) + 1; }
	  else {
	       qBegin = beginOffsetInQuery;
	       qEnd = endOffsetInQuery; }
	  int editDistance = 0;
	  for (int i=11; i<flds.size(); i++) {
	       if (strncmp (flds[i], "NM:i:", strlen ("NM:i:")) == 0) {
		    editDistance = atoi (flds[i]+5);
		    break; }
	  }
	  if ((isFirstRecord) || (strcmp (flds[0], querySeqNameHold) != 0) | (strcmp (flds[2], refSeqNameHold) != 0)) {
	       printf (">%s %s %d %d\n", flds[2], flds[0], len[std::string(flds[2])], currentOffset-1);
	       isFirstRecord = 0;
	       strcpy (refSeqNameHold, flds[2]);
	       strcpy (querySeqNameHold, flds[0]); }
	  printf ("%d %d %d %d %d %d %d\n", leftPos, endOffsetInReference, qBegin, qEnd, editDistance, editDistance, 0);
	  if (! noIndelLines)
	       printf ("%s", (char *) indelOutputSequence);
	  printf ("0\n");
     }

     return (0);
}


