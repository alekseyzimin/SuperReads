// The first argument is the file output by getUnitigTypeFromAsmFile.perl
// Each line of the output has two fields:
// (1) The unitig name, and (2) The unitig type, where the unitig type is
// 'U': U-unitig; 'S': surrogate; 'N': degenerate
// We are only interested in the surrogates
// The second arg is the utglen file
// The third arg is the frgutg file
// The fourth arg is the utgctg file
// The fifth arg is the frgctg file
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <charb.hpp>
#include <misc.hpp>

struct readPlacementStruct
{
     std::string readName;
     int beginOffset;
     int endOffset;
     char ori;
};
typedef std::vector<readPlacementStruct> readPlacementList;

bool readPlacementCompare (readPlacementStruct val1, readPlacementStruct val2);
FILE *Fopen (const char *fn, const char *mode);

int main (int argc, char **argv)
{
     FILE *infile;
     std::vector<char *> flds;
     charb line(100);
     std::set<std::string> keepUnitig;
     std::map<std::string, int> utgLens;
     std::map<std::string, readPlacementList> readPlacementsPerUnitig, unitigPlacementsPerContig,
	  readPlacementsPerContig;
     readPlacementStruct rPSRecord;
     std::map<std::string, std::string> surrogateUnitigToContigPlacement;

     infile = Fopen (argv[1], "r");
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  if (flds[1][0] == 'S')
	       keepUnitig.insert(std::string(flds[0]));
     }
     fclose (infile);

     infile = Fopen (argv[2], "r");
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  utgLens[std::string(flds[0])] = atoi (flds[1]); }
     fclose (infile);

     infile = Fopen (argv[3], "r");
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  if (keepUnitig.find (std::string(flds[1])) == keepUnitig.end())
	       continue;
	  rPSRecord.readName = std::string(flds[0]);
	  rPSRecord.beginOffset = atoi(flds[2]);
	  rPSRecord.endOffset = atoi(flds[3]);
	  rPSRecord.ori = flds[4][0];
	  readPlacementsPerUnitig[std::string(flds[1])].push_back (rPSRecord);
	  surrogateUnitigToContigPlacement[rPSRecord.readName] = std::string("S");
     }
     fclose (infile);

     infile = Fopen (argv[4], "r");
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  if (keepUnitig.find (std::string(flds[0])) == keepUnitig.end())
	       continue;
	  rPSRecord.readName = std::string(flds[0]); // Actually a unitig name here
	  rPSRecord.beginOffset = atoi(flds[2]);
	  rPSRecord.endOffset = atoi(flds[3]);
	  rPSRecord.ori = flds[4][0];
	  unitigPlacementsPerContig[std::string(flds[1])].push_back (rPSRecord);
     }
     fclose (infile);

     infile = Fopen (argv[5], "r");
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  rPSRecord.readName = std::string(flds[0]);
	  rPSRecord.beginOffset = atoi(flds[2]);
	  rPSRecord.endOffset = atoi(flds[3]);
	  rPSRecord.ori = flds[4][0];
	  readPlacementsPerContig[std::string(flds[1])].push_back (rPSRecord);
	  std::map<std::string, std::string>::iterator it =
	       surrogateUnitigToContigPlacement.find (rPSRecord.readName);
	  if (it != surrogateUnitigToContigPlacement.end())
	       surrogateUnitigToContigPlacement[rPSRecord.readName] =
		    std::string(flds[1]);
     }
     fclose (infile);

     // Now add the surrogates
     for (std::map<std::string, readPlacementList>::iterator it=unitigPlacementsPerContig.begin(); it!=unitigPlacementsPerContig.end(); ++it) {
	  std::string contigName = it->first;
	  // The following does one iteration per unitig in the contig
	  for (readPlacementList::iterator it2=(it->second).begin(); it2 != (it->second).end(); ++it2) {
	       readPlacementStruct unitigPlacementInContig = *it2;
	       std::string unitigName = unitigPlacementInContig.readName;
	       int unitigBegin = unitigPlacementInContig.beginOffset;
	       int unitigEnd = unitigPlacementInContig.endOffset;
	       char unitigOri = unitigPlacementInContig.ori;
	       // The following does one iteration per read in the unitig
	       for (readPlacementList::iterator it3=(readPlacementsPerUnitig[unitigName]).begin(); it3 != (readPlacementsPerUnitig[unitigName]).end(); ++it3) {
		    rPSRecord.readName = it3->readName;
		    if (surrogateUnitigToContigPlacement[rPSRecord.readName] == contigName)
			 continue;
#if 0
		    // This following is to not place if the read has already
		    // been placed in the frgctg file
		    if (surrogateUnitigToContigPlacement[rPSRecord.readName] != std::string("S"))
			 continue;
#endif
		    if (unitigOri == 'f') {
			 rPSRecord.beginOffset = unitigBegin + it3->beginOffset;
			 rPSRecord.endOffset = unitigBegin + it3->endOffset;
			 rPSRecord.ori = it3->ori; }
		    else {
			 rPSRecord.beginOffset = unitigEnd - it3->endOffset;
			 rPSRecord.endOffset = unitigEnd - it3->beginOffset;
			 if (it3->ori == 'f')
			      rPSRecord.ori = 'r';
			 else
			      rPSRecord.ori = 'f';
		    }
		    readPlacementsPerContig[contigName].push_back (rPSRecord);
	       }
	       // When we get here we've done all the reads in the unitig
	  }
	  // When we get here we've done all the unitigs for the contig
     }
     // When we get here we've added all the records needed for all contigs
	  
     std::vector<std::string> contigList;
     for (std::map<std::string, readPlacementList>::iterator it=readPlacementsPerContig.begin(); it!=readPlacementsPerContig.end(); ++it)
	  contigList.push_back (it->first);
     // Now we sort the contigs in numerical order
     sort (contigList.begin(), contigList.end());
//     for (int i=0; i<contigList.size(); i++) {
     for(auto contigName = contigList.cbegin(); contigName != contigList.cend(); ++contigName) {
//	  std::string contigName = contigList[i];
	  // Now sort the reads in placement order
	  std::sort (readPlacementsPerContig[*contigName].begin(), readPlacementsPerContig[*contigName].end(), readPlacementCompare);
	  
	  for (readPlacementList::iterator it=readPlacementsPerContig[*contigName].begin(); it!=readPlacementsPerContig[*contigName].end(); ++it) {
	       printf ("%s %s %d %d %c\n", it->readName.c_str(), contigName->c_str(), it->beginOffset, it->endOffset, it->ori);
	  }
     }
     
     return (0);
}

bool readPlacementCompare (readPlacementStruct val1, readPlacementStruct val2)
{
     if (val1.beginOffset != val2.beginOffset)
	  return (val1.beginOffset < val2.beginOffset);
     if (val1.ori != val2.ori)
	  return (val1.ori < val2.ori);
     if (val1.endOffset != val2.endOffset)
	  return (val1.endOffset < val2.endOffset);
     return (val1.readName < val2.readName);
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

