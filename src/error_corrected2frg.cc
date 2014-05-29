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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <charb.hpp>
#include <jellyfish/err.hpp>
#include <src/error_corrected2frg_cmdline.hpp>

namespace err = jellyfish::err;

inline void print_sequence(FILE *out, const char *seq, uint64_t len) {
  fprintf(out, "%s\n", seq);
  // uint64_t i;
  // for(i = 0; i < len; i += 70) {
  //   fprintf(out, "%.70s\n", seq + i);
  // }
  // if(i < len)
  //   fprintf(out, "%s\n", seq + i);
}

void print_quals(FILE *out, uint64_t len) {
  static const char *quals =
    "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";
  static const size_t qlen = strlen(quals);

  uint64_t i = 0;
  for( ; i + qlen < len; i += qlen) {
    fprintf(out, "%s", quals);
  }
  if(i < len)
    fprintf(out, "%.*s", (int)(len - i), quals);
  fprintf(out, "\n");
  
}

int main(int argc, char *argv[])
{
  FILE                        *FP;
  charb                        header;
  charb                        seq;
  std::vector<bool>            read_status;
  error_corrected2frg_cmdline  args(argc, argv);
  
  read_status.resize(args.number_of_reads_arg + 2);

  FP = fopen(args.fasta_file_arg,"r");
  if(!FP)
    err::die(err::msg() << "Can't open input file '" << args.fasta_file_arg << "': " << err::no);

  int first_char = fgetc(FP);
  if(first_char != '>')
    err::die("Invalid fasta file: no header found at beginning");

  printf("{VER\nver:2\n}\n"
         "{LIB\n" "act:A\n" "acc:%s\n"
         "ori:I\n" "mea:%.3g\n" "std:%.3g\n"
         "src:\n.\n" "nft:1\n" "fea:\n"
         "doNotOverlapTrim=1\n.\n"
         "}\n",
         args.lib_id_arg, args.mean_arg, args.stdev_arg);

   // First pass through reads: mark reads with "good" length
  while(fgets(header, FP)) {
     // Find read number
    const char *first_word = strtok(header, " \n");
    if(!first_word)
      err::die("Invalid empty header");
    uint64_t n = atoll(first_word+2);
 
    // Compute length of sequence
    uint64_t seq_len = 0;
    first_char = fgetc(FP);
    while(first_char != EOF && first_char != '>') {
      fgets(seq, FP);
      seq.chomp();
       seq_len += strlen(seq) + 1;
      first_char = fgetc(FP);
     }
     if(seq_len >= args.length_min_arg)
      read_status[n] = true;
   }

  // Second pass through reads
  rewind(FP);
  first_char = fgetc(FP);
  while(fgets(header, FP)) {
    const char *first_word = strtok(header, " \n");
    if(first_word == NULL)
      err::die("Invalid empty header in input fasta file");
    uint64_t n = atoll(first_word+2);

    first_char = fgetc(FP);
    seq.clear();
    while(first_char != EOF && first_char != '>') {
      ungetc(first_char, FP);
      fgets_append(seq, FP);
      seq.chomp();
      first_char = fgetc(FP);
    }
    
    // Check that we have a pair of reads
    if(!read_status[n] || !read_status[n + (n%2 == 0 ? 1 : -1)])
      continue;
    
    uint64_t seq_len = strlen(seq);
    printf("{FRG\n" "act:A\n" "acc:%s\n"
           "rnd:1\n" "sta:G\n" "lib:%s\n"
           "pla:0\n" "loc:0\n" "src:\n.\n"
           "seq:\n",
           (char*)header, args.lib_id_arg);
    print_sequence(stdout, seq, seq_len);
    printf(".\nqlt:\n");
    print_quals(stdout, seq_len);
    printf(".\n" "hps:\n.\n" "clr:0,%ld\n"
           "}\n", seq_len);
  }

  // Print linkage messages
  for(uint64_t i = 0; i <= args.number_of_reads_arg; i += 2) {
    if(read_status[i] && read_status[i+1]) {
      printf("{LKG\nact:A\nfrg:%s%ld\nfrg:%s%ld\n}\n", 
             args.lib_id_arg, i, args.lib_id_arg, i+1);
    }
  }

  fclose(FP);
  return(0);
}

