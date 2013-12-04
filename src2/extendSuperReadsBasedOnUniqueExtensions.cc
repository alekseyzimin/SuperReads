#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <misc.hpp>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <src2/extendSuperReadsBasedOnUniqueExtensions_cmdline.hpp>

std::string orientSuperRead (std::string inputString);
FILE *Fopen (const char *fn, const char *mode);

struct stringAndLenStruct {
     std::string str;
     int len; };

bool stringAndLenCompare (struct stringAndLenStruct i, struct stringAndLenStruct j)
{
     if (i.len != j.len)
	  return (i.len > j.len);
     return (i.str < j.str);
}

int main (int argc, char **argv)
{
     charb dir;
     cmdline_parse args;
     args.parse (argc, argv);
     strcpy (dir, args.dir_arg);
     int overlap = args.mer_arg;
     charb fname(100);
     charb line(100);
     std::vector<char *> flds, flds2;
     int maxUnitigNum;
     sprintf (fname, "%s/maxKUnitigNumber.txt", (char *) dir);
     FILE *infile = Fopen (fname, "r");
     int numItems = fscanf (infile, "%d\n", &maxUnitigNum);
     if (numItems != 1)
	  exit (1);
     fclose (infile);
     int *afterUni = (int *) calloc (maxUnitigNum+1, sizeof(int));
     char *afterUniOri = (char *) calloc (maxUnitigNum+1, sizeof(char));
     unsigned char *afterUniBad = (unsigned char *) calloc (maxUnitigNum+1, sizeof (unsigned char ));
     int *beforeUni = (int *) calloc (maxUnitigNum+1, sizeof(int));
     char *beforeUniOri = (char *) calloc (maxUnitigNum+1, sizeof(char));
     unsigned char *beforeUniBad = (unsigned char *) calloc (maxUnitigNum+1, sizeof (unsigned char ));
     unsigned int *alreadySeen = (unsigned int *) calloc (maxUnitigNum+1, sizeof (unsigned int));
     for (int i=0; i<=maxUnitigNum; ++i)
	  afterUni[i] = beforeUni[i] = -1;
     sprintf (fname, "%s/extendSuperReadsForUniqueKmerNeighbors.output.txt", (char *) dir);
     infile = Fopen (fname, "r");
     while (fgets (line, 100, infile)) {
	  int numFlds = getFldsFromLine (line, flds);
	  int uni = atoi (flds[1]);
	  char ori = flds[2][strlen(flds[2])-1];
	  flds[2][strlen(flds[2])-1] = 0;
	  int otherUni = atoi (flds[2]);
	  unsigned char isBad;
	  if (numFlds == 3)
	       isBad = 0;
	  else
	       isBad = 1;
	  if (strcmp (flds[0], "BEF") == 0) {
	       beforeUni[uni] = otherUni;
	       beforeUniOri[uni] = ori;
	       beforeUniBad[uni] = isBad; }
	  else {
	       afterUni[uni] = otherUni;
	       afterUniOri[uni] = ori;
	       afterUniBad[uni] = isBad; }
     }
     fclose (infile);

     sprintf (fname, "%s/kUnitigLengths.txt", (char *) dir);
     int *unitigLengths = (int *) calloc (maxUnitigNum+1, sizeof(int));
     infile = Fopen (fname, "r");
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  int uni = atoi (flds[0]);
	  unitigLengths[uni] = atoi (flds[1]); }
     fclose (infile);

//     sprintf (fname, "%s/reduce.tmp", (char *) dir);
     std::vector<std::string> listOfInputSuperReads;
     std::vector<int> inputSuperReadLengths;
     std::set<std::string> superReadAlreadyIn;
     std::vector<stringAndLenStruct> listOfNewSuperReads;
     sprintf (fname, "%s/sr_sizes.tmp.hold", (char *) dir);
     infile = Fopen (fname, "r");
     unsigned int currentLineValue = 0;
     int fldNumToUse = 0;
     int superReadLength;
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  std::string superReadName = std::string(flds[fldNumToUse]);
	  superReadLength = atoi (flds[1]);
          listOfInputSuperReads.push_back (superReadName);
          inputSuperReadLengths.push_back (superReadLength);
	  superReadAlreadyIn.insert (superReadName);
	  int numFlds = getFldsFromLine (flds[fldNumToUse], flds2, "_");
	  ++currentLineValue;
	  for (int i=0; i<numFlds; ++i)
	       alreadySeen[atoi(flds2[i])] = currentLineValue;
	  std::vector<int> preUnis, postUnis;
	  std::vector<char> preOris, postOris;
	  char firstOri = flds2[0][strlen(flds2[0])-1], lastOri;
	  flds2[0][strlen(flds2[0])-1] = 0;
	  int firstUni = atoi (flds2[0]), lastUni;
	  if (numFlds == 1) {
	       lastUni = firstUni;
	       lastOri = firstOri; }
	  else {
	       lastOri = flds2[numFlds-1][strlen(flds2[numFlds-1])-1];
	       flds2[numFlds-1][strlen(flds2[numFlds-1])-1] = 0;
	       lastUni = atoi (flds2[numFlds-1]); }
	  std::vector<int> beforeUnitigs, afterUnitigs;
	  std::vector<char> beforeOris, afterOris;
	  while (1) {
	       char tempOri;
	       int tempUni;
//	       printf ("In part 1, superRead = %s\n", superReadName.c_str());
	       if (firstOri == 'F') {
		    tempUni = beforeUni[firstUni];
		    if (beforeUniBad[firstUni])
			 break;
		    tempOri = beforeUniOri[firstUni];
	       }
	       else {
		    tempUni = afterUni[firstUni];
		    if (afterUniBad[firstUni])
			 break;
		    if (afterUniOri[firstUni] == 'F')
			 tempOri = 'R';
		    else
			 tempOri = 'F';
	       }
	       if (tempUni < 0)
		    break;
	       if (alreadySeen[tempUni] == currentLineValue)
		    break;
	       alreadySeen[tempUni] = currentLineValue;
	       beforeUnitigs.push_back (tempUni);
	       beforeOris.push_back (tempOri);
	       firstUni = tempUni;
	       firstOri = tempOri;
	  }
	  while (1) {
	       char tempOri;
	       int tempUni;
//	       printf ("In part 2, superRead = %s\n", superReadName.c_str());
	       if (lastOri == 'F') {
		    tempUni = afterUni[lastUni];
		    if (afterUniBad[lastUni])
			 break;
		    tempOri = afterUniOri[lastUni];
	       }
	       else {
		    tempUni = beforeUni[lastUni];
		    if (beforeUniBad[lastUni])
			 break;
		    if (beforeUniOri[lastUni] == 'F')
			 tempOri = 'R';
		    else
			 tempOri = 'F';
	       }
	       if (tempUni < 0)
		    break;
	       if (alreadySeen[tempUni] == currentLineValue)
		    break;
	       alreadySeen[tempUni] = currentLineValue;
	       afterUnitigs.push_back (tempUni);
	       afterOris.push_back (tempOri);
	       lastUni = tempUni;
	       lastOri = tempOri;
	  }
	  // If nothing extending the ends then nothing to do
	  if (beforeUnitigs.size() + afterUnitigs.size() == 0) {
//	       printf ("%s\n", superReadName.c_str());
	       continue; }
	  std::reverse (beforeUnitigs.begin(), beforeUnitigs.end());
	  std::reverse (beforeOris.begin(), beforeOris.end());
	  int newStringLength = superReadLength;
//	  printf ("%s ", superReadName.c_str());
	  std::string newSuperReadString;
	  for (unsigned int i=0; i<beforeUnitigs.size(); ++i) {
	       charb uniNum(30);
	       sprintf (uniNum, "%d", beforeUnitigs[i]);
	       newSuperReadString += (char *) uniNum;
	       newSuperReadString += beforeOris[i];
	       newSuperReadString += "_";
	       newStringLength += (unitigLengths[beforeUnitigs[i]] - overlap);
//	       printf ("%d%c ", beforeUnitigs[i], beforeOris[i]);
	  }
//	  printf ("%s", superReadName.c_str());
	  newSuperReadString += superReadName;
	  for (unsigned int i=0; i<afterUnitigs.size(); ++i) {
	       charb uniNum(30);
	       sprintf (uniNum, "%d", afterUnitigs[i]);
	       newSuperReadString += "_";
	       newSuperReadString += (char *) uniNum;
	       newSuperReadString += afterOris[i];
	       newStringLength += (unitigLengths[afterUnitigs[i]] - overlap);
//	       printf (" %d%c", afterUnitigs[i], afterOris[i]);
	  }
//	  printf ("\n");
	  if (superReadAlreadyIn.find (newSuperReadString) != superReadAlreadyIn.end())
	       continue; // Already found
	  newSuperReadString = orientSuperRead (newSuperReadString);
	  if (superReadAlreadyIn.find (newSuperReadString) != superReadAlreadyIn.end())
	       continue; // Already found
	  // Add it to the new super-reads along with its length
	  struct stringAndLenStruct stdAndLenRec;
	  stdAndLenRec.str = newSuperReadString;
	  stdAndLenRec.len = newStringLength;
	  listOfNewSuperReads.push_back (stdAndLenRec);
	  // Then mark it as added
	  superReadAlreadyIn.insert (newSuperReadString);
//	  printf ("\nnewSuperRead = %s\n", newSuperReadString.c_str());
     }
     fclose (infile);

     sort (listOfNewSuperReads.begin(), listOfNewSuperReads.end(), stringAndLenCompare);
     fprintf (stderr, "listOfNewSuperReads.size() = %d\n", (int) listOfNewSuperReads.size());

#if 0
     printf ("New superReads:\n");
     for (unsigned int i=0; i<listOfNewSuperReads.size(); ++i) {
	  printf ("%s %d\n", listOfNewSuperReads[i].str.c_str(), listOfNewSuperReads[i].len); }
#endif
     int maxLen = inputSuperReadLengths[0];
     if (listOfNewSuperReads[0].len > maxLen)
	  maxLen = listOfNewSuperReads[0].len;
     unsigned int oldIndex=0, newIndex=0;
//     printf ("All superReads together:\n");
     sprintf (fname, "%s/superReadNames.addedByExtension.txt", (char *) dir);
     FILE *outfile = Fopen (fname, "w");
     for (int len=maxLen; len>0; --len) {
	  while ((oldIndex<inputSuperReadLengths.size()) && (inputSuperReadLengths[oldIndex] == len)) {
	       printf ("%s %d\n", listOfInputSuperReads[oldIndex].c_str(), inputSuperReadLengths[oldIndex]);
	       ++oldIndex; }
	  while ((newIndex<listOfNewSuperReads.size()) && (listOfNewSuperReads[newIndex].len == len)) {
	       printf ("%s %d\n", listOfNewSuperReads[newIndex].str.c_str(), listOfNewSuperReads[newIndex].len);
	       fprintf (outfile, "%s\n", listOfNewSuperReads[newIndex].str.c_str());
	       ++newIndex; }
     }
     fclose (outfile);

     return (0);
}

std::string orientSuperRead (std::string inputString)
{
     int numFlds;
     std::vector<char *>flds;
     charb localLine(100);
     std::vector<int> uniNums;
     std::vector<char> oris;
     strcpy (localLine, inputString.c_str());
     numFlds = getFldsFromLine (localLine, flds, "_");
//     printf ("inputString = %s\n", inputString.c_str());
     fflush (stdout);
     for (int i=0; i<numFlds; ++i) {
//	  printf ("flds[%d] = %s\n", i, flds[i]); fflush (stdout);
	  uniNums.push_back (atoi (flds[i]));
	  oris.push_back (flds[i][strlen(flds[i])-1]);
//	  printf ("%s ", flds[i]); 
     }
     unsigned char stringAction = 2;
     for (int i=0, j=numFlds-1; i<=j; ++i, --j) {
	  if (uniNums[i] != uniNums[j]) {
	       if (uniNums[i] < uniNums[j])
		    stringAction = 0;
	       else
		    stringAction = 1;
	       break; }
	  if (oris[i] == oris[j]) {
	       if (oris[i] == 'F')
		    stringAction = 0;
	       else
		    stringAction = 1;
	       break; }
	  if (stringAction != 2)
	       break;
     }
//     printf ("\n");
     
     if (stringAction == 0)
	  return (inputString);
     else {
	  std::string outputString;
	  sprintf (localLine, "%d%c", uniNums[numFlds-1], (oris[numFlds-1] == 'F') ? 'R' : 'F');
	  
	  for (int i=numFlds-2; i>=0; --i)
	       sprintf_append (localLine, "_%d%c", uniNums[i], (oris[i] == 'F') ? 'R' : 'F');
	  outputString = std::string (localLine);
	  return (outputString);
     }
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

