// For now make the jumping-read fasta file fixed and the placement file hardcoded
// Later the reads will be gotten from stdin and the placement file via arg
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <charb.hpp>
#include <misc.hpp>
#include <src2/putReadsIntoGroupsBasedOnSuperReads_cmdline.hpp>
int main (int argc, char **argv)
{
     charb line(100), cmd(100);
     FILE *infile;
     std::vector<std::vector<std::string> > superReadToReadList;
     cmdline_parse args;
     args.parse (argc, argv);
     sprintf (cmd, "tail -2 %s | head -1", args.super_read_sequence_file_arg);
     infile = popen ((char *)cmd, "r");
//     infile = popen ("tail -2 /genome3/raid/tri/reduceReadsIntoCelera/assembly/work2/superReadSequences.fasta | head -1", "r");
     fgets (line, 100, infile);
     int lastSuperReadNum = atoi (line+1);
     fclose (infile);
     std::vector<std::string> emptyStringVector;
     for (int i=0; i<=lastSuperReadNum; ++i)
	  superReadToReadList.push_back (emptyStringVector);
//     infile = fopen ("/genome3/raid/tri/reduceReadsIntoCelera/assembly/sj.cor.clean.fa", "r");
     infile = stdin;
     charb rname(100);
     std::map<std::string, std::string> readNamesToReadSequences;
     while (fgets (line, 100, infile)) {
	  sscanf ((char *) line, ">%s", (char *) rname);
	  for (int i=1; i<4; ++i)
	       fgets_append (line, infile);
	  readNamesToReadSequences[std::string ((char *) rname)] = std::string(line);
     }
//     fclose (infile);
     std::vector<char *> flds;
     infile = fopen (args.read_placements_file_arg, "r");
//     infile = fopen ("/genome3/raid/tri/reduceReadsIntoCelera/assembly/work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.reduced.txt", "r");
     std::set<std::string> wasFound;
     std::map<std::string, int> superReadForReadName;
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  if (flds[0][strlen(flds[0])-1] % 2 == 1)
	       --flds[0][strlen(flds[0])-1];
	  std::string readName = std::string(flds[0]);
	  int superRead = atoi (flds[1]);
	  if (wasFound.find (readName) == wasFound.end()) {
	       superReadForReadName[readName] = superRead;
	       wasFound.insert (readName);
	       continue; }
	  // If we get here it's the second of the pair
	  int superRead2 = atoi (flds[1]);
	  if (superRead2 > superReadForReadName[readName])
	       superRead = superReadForReadName[readName];
	  superReadToReadList[superRead].push_back (readName); }
     fclose (infile);
     
     for (int i=0; i<=lastSuperReadNum; ++i)
	  for (unsigned int j=0; j<superReadToReadList[i].size(); ++j)
	       std::cout << readNamesToReadSequences[superReadToReadList[i][j]];
     
     return (0);
}

