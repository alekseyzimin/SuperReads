#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
using namespace std;

char *flds[10000];
int getFldsFromLine (char *cptrHold);

int main (int argc, char **argv)
{
     map<string,int> superReadToCounts;
     map<string,int>::iterator it;
     string str;
     int isFirstLine, numFlds;
     char line[1000000], line2[1000000], *cptr;
     char *prefixHold;
     long long readNumHold, readNum;
     char *superReadHold;

     isFirstLine = 1;
     cptr = line;
     fgets (cptr, 1000000, stdin);
     numFlds = getFldsFromLine (cptr);
     str = string (flds[1]);
     superReadToCounts.insert (superReadToCounts.begin(), pair<string,int>(str, 1));
     superReadHold = flds[1];
     prefixHold = flds[0]; // Just use the first 2 chars
     readNumHold = atoll (cptr+2);
     cptr = line2;
     while (fgets (cptr, 1000000, stdin)) {
	  numFlds = getFldsFromLine (cptr);
	  readNum = atoll(cptr+2);
	  // The next excludes counting the second read of a mate pair from
	  // the counts for the super-read. We want super-reads which have
	  // reads coming from at least 2 different inserts. Our read counts
	  // for the inserts will (probably) be off.
#if 1
	  if ((strcmp(flds[1], superReadHold) == 0) &&
	      (readNumHold+1 == readNum) &&
	      (readNum % 2 == 1) &&
	      (strncmp(prefixHold, cptr, 2) == 0))
	       continue;
#endif
	  str = string (flds[1]);
	  it = superReadToCounts.find(str);
	  if (it == superReadToCounts.end()) {
	       it = superReadToCounts.begin();
	       superReadToCounts.insert (it, pair<string,int>(str, 1));
	  }
	  else {
	       ++(*it).second;
	  }
	  superReadHold = flds[1];
	  prefixHold = flds[0];
	  readNumHold = readNum;
	  if (cptr == line)
	       cptr = line2;
	  else
	       cptr = line;
     }
     for (it=superReadToCounts.begin(); it != superReadToCounts.end(); it++)
	  cout << (*it).second << " " << (*it).first << endl;

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

