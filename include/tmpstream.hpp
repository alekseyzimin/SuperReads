#ifndef __TMP_FILE_H_
#define __TMP_FILE_H_

#include <ext/stdio_filebuf.h>
#include <stdio.h>
#include <stdexcept>

/** Temporary file stream. Wrapper around FILE* tmpfile(). It seems
 * that stdio_filebuf creates a memory leak by not closing fclose on
 * the file descriptor. Oh well!
 */

template<typename _CharT, typename _Traits = std::char_traits<_CharT> >
class basic_tmpstream : public std::iostream {
  typedef __gnu_cxx::stdio_filebuf<_CharT> stdbuf;
public:
  basic_tmpstream() : std::iostream(open_tmp()) { }
  virtual ~basic_tmpstream() { close(); delete rdbuf(); }
  bool is_open() { return ((stdbuf*)rdbuf())->is_open(); }
  void close() { ((stdbuf*)rdbuf())->close(); }
private:
  static stdbuf * open_tmp() {
    FILE *f = tmpfile();
    if(!f)
      throw std::runtime_error("Failed to create a temporary file.");
    return new stdbuf(f, std::ios::in|std::ios::out);
  }
};
typedef basic_tmpstream<char> tmpstream;


#endif /* __TMP_FILE_H_ */
