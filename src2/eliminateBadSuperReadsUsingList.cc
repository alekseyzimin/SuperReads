/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/* This program expects 2 args:
   1) The name of the file containing the read placements in the super-reads
   2) The file containing passing super read names produced when creating the super-read sequences
   The file containing read placements for reads in the good super-reads
      will come out on stdout.
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <set>
#include <vector>
#include <charb.hpp>
#include <map>
#include <iostream>
#include <misc.hpp>
using namespace std;

struct oldSuperReadPlacementStruct
{
     string newSuperReadName;
     int newOffsetToStart;
     char newOri;
};

charb line(1000000);
vector<char *> flds;
map<string, struct oldSuperReadPlacementStruct> oldSuperReadToNewSuperReadMap;

FILE *Fopen (const char *fn, const char *mode);

int main (int argc, char **argv)
{
     set<string> isGoodSuperRead;
     FILE *infile;
     char *goodFilename=NULL;
     char *readPlacementFilename = NULL;
     char *superReadReductionFile = NULL;
     char *cptr;
     for (int i=1; i<argc; i++) {
	  if (strcmp(argv[i], "--reduce-file") == 0) {
	       ++i;
	       superReadReductionFile = argv[i];
	       continue; }
          if (strcmp(argv[i], "--good-super-reads-file") == 0) {
               ++i;
               goodFilename = argv[i];
               continue; }
          if (strcmp(argv[i], "--read-placement-file") == 0) {
               ++i;
               readPlacementFilename = argv[i];
               continue; }
          //if we get here, then there is a problem with the arguments
               fprintf (stderr, "Unrecognized flag '%s'. Bye!\n", argv[i]);
               exit (1); 
     }

     if(readPlacementFilename == NULL){
	fprintf (stderr, "No input read placement file specified. Bye!\n");
        exit (1); 
     }
       
     if (superReadReductionFile != NULL) {
	  infile = Fopen (superReadReductionFile, "r");
	  while (fgets (line, 1000000, infile)) {
	       oldSuperReadPlacementStruct osrps;
	       getFldsFromLine (line, flds);
	       osrps.newSuperReadName = string(flds[1]);
	       osrps.newOri = *flds[2];
	       osrps.newOffsetToStart = atoi(flds[3]);
	       oldSuperReadToNewSuperReadMap[string(flds[0])] = osrps;
	  }
	  fclose (infile);
     }
     
     if (goodFilename != NULL) {
     infile = Fopen (goodFilename, "r");
     while (fgets (line, 1000000, infile)) {
	  if (! isdigit(line[0]))
	       continue;
	  cptr = line;
	  while (! isspace(*cptr)) ++cptr;
	  *cptr = 0;
	  isGoodSuperRead.insert (string(line));
     }
     fclose (infile);
     }

     infile = Fopen (readPlacementFilename, "r");
     while (fgets (line, 1000000, infile)) {
	  getFldsFromLine (line, flds);
	  string superRead = string(flds[1]);
          if(goodFilename != NULL){
	  if (isGoodSuperRead.find (superRead) == isGoodSuperRead.end())
	       continue;
          }
	  if (superReadReductionFile == NULL) {
	       for (int i=0; i<4; i++) {
		    fputs (flds[i], stdout);
		    if (i<3)
			 fputc (' ', stdout);
		    else
			 fputc ('\n', stdout); }
	       continue; }
	       
	  // We only get here if there is a superReadReductionFile
	  auto it = oldSuperReadToNewSuperReadMap.find(superRead);
	  if (it == oldSuperReadToNewSuperReadMap.end()) {
	       for (int i=0; i<4; i++) {
		    fputs (flds[i], stdout);
		    if (i<3)
			 fputc (' ', stdout);
		    else
			 fputc ('\n', stdout); }
	       continue; }
	  int offset = atoi(flds[2]);
	  char ori = *flds[3];
	  if (it->second.newOri == 'F') {
	       offset += (it->second).newOffsetToStart; }
	  else {
	       ori = (ori == 'F') ? 'R' : 'F';
	       offset = (it->second).newOffsetToStart - offset;
	  }
	  cout << flds[0] << " " << (it->second).newSuperReadName << " " << offset << " " << ori << '\n';
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

