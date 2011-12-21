#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <exp_buffer.hpp>
#include <charb.hpp>
#include <err.hpp>
#include <misc.hpp> // for getFldsFromLine
#include <src/sr_names.hpp>
#include <src2/getSuperReadInsertCountsFromReadPlacementFile.hpp>

using namespace std;

// The following feature is broken! The char * for prefixHold and
// superReadHold need to be changed to a charb as well (i.e., to copy
// the string from one iteration to the next).
/*  #define EXCLUDE_MATE 1 */

struct str_comp {
  bool operator()(const char *s1, const char *s2) {
    return strcmp(s1, s2) < 0;
  }
};

typedef const char*(*encode_fn)(const char* str);
typedef const char*(*decode_fn)(const char* str);

const char* identity(const char* str) { return str; }
const char* dup(const char* str) { return strdup(str); }
const char* fib_decode(const char* str) {
  static charb buffer;
  sr_name::to_str(str, buffer);
  return buffer;
}

int main (int argc, char **argv)
{
  getSuperReadInsertCountsFromReadPlacementFile args(argc, argv);

  typedef map<const char *, int, str_comp> read_count_map;
  read_count_map    superReadToCounts;
  charb*            lines = new charb[2];
  charb*            cptr;
  int               line  = 1;
#if EXCLUDE_MATE
  char*             prefixHold;
  long long         readNumHold, readNum;
  char*             superReadHold;
#endif
  ExpBuffer<char *> flds;
  std::ifstream input(args.input_arg);
  if(!input.good())
    die << "Can't open input file '" << args.input_arg << "'" << err::no;
  std::ofstream output(args.output_arg);
  if(!output.good())
    die << "Can't open output file '" << args.output_arg << "'" << err::no;

  encode_fn encode;
  decode_fn decode;

  if(args.fib_flag) {
    encode = sr_name::from_str;
    decode = fib_decode;
  } else {
    encode = dup;
    decode = identity;
  }

  cptr = &lines[(line = !line)];
  if(!getline (input, *cptr))
    die << "Invalid input: empty file";
  getFldsFromLine (*cptr, flds);
  const char *nstr = encode(flds[1]);
  ++superReadToCounts[nstr];
#if EXCLUDE_MATE
  superReadHold = flds[1];
  prefixHold    = flds[0];      // Just use the first 2 chars
  readNumHold   = atoll (*cptr+2);
#endif
  cptr          = &lines[(line = !line)];
  while (getline (input, *cptr)) {
    getFldsFromLine (*cptr, flds);
#if EXCLUDE_MATE
    // The next excludes counting the second read of a mate pair from
    // the counts for the super-read. We want super-reads which have
    // reads coming from at least 2 different inserts. Our read counts
    // for the inserts will (probably) be off.
    readNum = atoll(*cptr+2);
    if ((strcmp(flds[1], superReadHold) == 0) &&
        (readNumHold+1 == readNum) &&
        (readNum % 2 == 1) &&
        (strncmp(prefixHold, *cptr, 2) == 0))
      continue;
#endif
    nstr = encode(flds[1]);
    ++superReadToCounts[nstr];
#if EXCLUDE_MATE
    superReadHold = flds[1];
    prefixHold    = flds[0];
    readNumHold   = readNum;
#endif
    cptr          = &lines[(line = !line)];
  }
  for (read_count_map::const_iterator it = superReadToCounts.begin(); 
       it != superReadToCounts.end(); ++it)
    output << it->second << " " << decode(it->first) << "\n";

  return (0);
}
