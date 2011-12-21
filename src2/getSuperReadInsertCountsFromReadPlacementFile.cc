#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <err.hpp>
#include <misc.hpp> // for getFldsFromLine

using namespace std;

struct str_comp {
  bool operator()(const char *s1, const char *s2) {
    return strcmp(s1, s2) < 0;
  }
};

int main (int argc, char **argv)
{
  typedef map<const char *, int, str_comp> read_count_map;
  read_count_map    superReadToCounts;
  charb*            lines = new charb[2];
  charb*            cptr;
  int               line  = 1;
  char*             prefixHold;
  long long         readNumHold, readNum;
  char*             superReadHold;
  ExpBuffer<char *> flds;

  cptr = &lines[(line = !line)];
  if(!fgets (*cptr, stdin))
    die << "Invalid input: empty file";
  getFldsFromLine (*cptr, flds);
  const char *nstr = strdup(flds[1]);
  ++superReadToCounts[nstr];
  superReadHold = flds[1];
  prefixHold    = flds[0];      // Just use the first 2 chars
  readNumHold   = atoll (*cptr+2);
  cptr          = &lines[(line = !line)];
  while (fgets (*cptr, stdin)) {
    getFldsFromLine (*cptr, flds);
    readNum = atoll(*cptr+2);
    // The next excludes counting the second read of a mate pair from
    // the counts for the super-read. We want super-reads which have
    // reads coming from at least 2 different inserts. Our read counts
    // for the inserts will (probably) be off.
#if 1
    if ((strcmp(flds[1], superReadHold) == 0) &&
        (readNumHold+1 == readNum) &&
        (readNum % 2 == 1) &&
        (strncmp(prefixHold, *cptr, 2) == 0))
      continue;
#endif
    nstr = strdup(flds[1]);
    ++superReadToCounts[nstr];
    superReadHold = flds[1];
    prefixHold    = flds[0];
    readNumHold   = readNum;
    cptr          = &lines[(line = !line)];
  }
  for (read_count_map::const_iterator it = superReadToCounts.begin(); 
       it != superReadToCounts.end(); ++it)
    cout << it->second << " " << it->first << "\n";

  return (0);
}
