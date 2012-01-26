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
#include <src2/getSuperReadInsertCountsFromReadPlacementFileTwoPasses.hpp>
#include <src/bloom_counter2.hpp>

using namespace std;

#define EXCLUDE_MATE 1

struct str_comp {
  bool operator()(const char *s1, const char *s2) {
    return strcmp(s1, s2) < 0;
  }
};

typedef const char*(*coding_fn)(const char* str);

// Function pointer to encode (if so desired) the names.
const char* str_identity(const char* str) { return str; }
const char* str_dup(const char* str) { return strdup(str); }
const char* fib_decode(const char* str) {
  static charb buffer;
  sr_name::to_str(str, buffer);
  return buffer;
}

// Parse the input input stream is. Delegate to s the decision of
// where to store the read name
template<typename store>
void parse_store_input(std::istream &is, store &s) {
  const char*       superReadHold = "";
  const char*       prefixHold    = "";
  long long         readNumHold   = 0;
  
  ExpBuffer<char *> flds;
  charb*            lines = new charb[2];
  int               line  = 0;
  charb*            cptr  = &lines[line];

  while(getline(is, *cptr)) {
    getFldsFromLine(*cptr, flds);
    long long readNum = atoll(*cptr+2);
    if(strcmp(flds[1], superReadHold) == 0 &&
       readNumHold+1 == readNum &&
       readNum % 2 == 1 &&
       strncmp(prefixHold, *cptr, 2) == 0)
      continue; // Exclude second read from a mate-pair
    s(flds[1]);
    superReadHold = flds[1];
    prefixHold    = flds[0];
    readNumHold   = readNum;
    cptr = &lines[(line = !line)];
  }
}

// Store in a bloom filter
struct bloom_store {
  bloom_counter2<const char*> bc;
  bloom_store(size_t n) : bc(0.01, n) { }
  void operator()(const char *s) { bc.insert(s); }
};

// Store in a map provided that the entry has been seen more than once
struct map_store {
  bloom_counter2<const char*>& bc;
  typedef std::map<const char*, int, str_comp> map_type;
  typedef map_type::const_iterator iterator;
  map_type map;
  coding_fn encode;
  map_store(bloom_counter2<const char*>& bc_, coding_fn encode_) :
    bc(bc_), map(), encode(encode_) { }
  void operator()(const char *s) {
    if(bc.check(s) > (unsigned int)1)
      ++map[encode(s)];
  }
};

int main (int argc, char **argv)
{
  typedef getSuperReadInsertCountsFromReadPlacementFileTwoPasses arg_parse;
  arg_parse args(argc, argv);

  std::ofstream output(args.output_arg);
  if(!output.good())
    die << "Can't open output file '" << args.output_arg << "'" << err::no;

  // Parse input into bloom counter
  bloom_store bs(args.number_reads_arg);
  for(arg_parse::input_arg_const_it file = args.input_arg.begin(); file != args.input_arg.end(); ++file) {
    std::ifstream input(*file);
    parse_store_input(input, bs);
  }

  // Parse input into map, if count > 1
  coding_fn encode = args.fib_flag ? (coding_fn)sr_name::from_str : str_dup;
  map_store ms(bs.bc, encode);
  for(arg_parse::input_arg_const_it file = args.input_arg.begin(); file != args.input_arg.end(); ++file) {
    std::ifstream input(*file);
    parse_store_input(input, ms);
  }

  // Output result
  coding_fn decode = args.fib_flag ? fib_decode : str_identity;
  for (map_store::iterator it = ms.map.begin(); it != ms.map.end(); ++it)
    if(it->second > 1)
      output << it->second << " " << decode(it->first) << "\n";
  output.close();

  return (0);
}
