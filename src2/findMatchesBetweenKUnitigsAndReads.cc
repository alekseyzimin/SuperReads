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


/* Arguments are as follows:
   1) The name of the database output by jellyfish (using raw data)
   2) The name of the file with the k-unitigs
   3) The name of the file which reports how many k-unitigs there are
   4) The name of the file containing the reads

   4/12/11: The output has now been reduced; to get the longer (former output),
   use the '-l' flag.
   5/11/11: Made parallel
   5/11/11: Use the '-p' flag to give an output prefix.
   5/13/11: Use the '-t' flag to specify how many threads should run.
*/

#ifndef typeof
#define typeof __typeof__
#endif

#include <multi_thread_skip_list_map.hpp>
#include <jellyfish/err.hpp>
#include <src/mer_dna.hpp>
#include <src/read_parser.hpp>

#include <jellyfish/mer_counting.hpp>
#include <jellyfish/mapped_file.hpp>
#include <jellyfish/invertible_hash_array.hpp>
#include <jellyfish/allocators_mmap.hpp>
#include <stdint.h>
#include <cstdio>

#include <iostream>
#include <fstream>

#include <src/diskBasedUnitigger.h>
#define KUNITIG_FILE "/genome8/raid/tri/kUnitigStudy/arg_ant/afterAlekseyAndMikeRedundentKill/guillaumeKUnitigsAtLeast32bases_all.fasta"
#define READ_DATA_FILE "/genome8/raid/tri/testDirForReadPlacementRoutines/brucellaData/brucella.pass5reads.fasta"

#include <charb.hpp>
#include <gzip_stream.hpp>
#include <jflib/multiplexed_io.hpp>

#include <src2/findMatchesBetweenKUnitigsAndReads_cmdline.hpp>

struct kMerUnitigInfoStruct {
  uint32_t kUnitigNumber;
  uint32_t kUnitigOffset:31;
  uint32_t kmerOriInKunitig:1;
};

/* Our types.
 */
//typedef jellyfish::invertible_hash::array<uint64_t,atomic::gcc<uint64_t>,allocators::mmap> inv_hash_storage_t;
typedef inv_hash_storage_t::iterator iterator_t;
class header_output;
class mer_unitig_info_base {
public:
  virtual const kMerUnitigInfoStruct* find(const mer_dna& mer) = 0;
};

void getMatchesForRead (const char *readBasesBegin, const char *readBasesEnd, 
                        const char *readName, mer_unitig_info_base* mer_unitig_info,
                        header_output& out);


class jf_kmer_unitig_info : public mer_unitig_info_base {
  typedef ExpBuffer<kMerUnitigInfoStruct, remaper<kMerUnitigInfoStruct>> kMerUnitigInfo_type;
  mapped_file          dbf_;
  raw_inv_hash_query_t qhash_;
  inv_hash_storage_t*  hash_;
  kMerUnitigInfo_type  kMerUnitigInfo_;

public:
  jf_kmer_unitig_info(const char* hash_path) :
    dbf_(hash_path), qhash_(dbf_), hash_(qhash_.get_ary()), kMerUnitigInfo_(hash_->get_size())
  { 
    kMerUnitigInfo_.touch_all();
    dbf_.random().load();
  }

  unsigned int mer_len() { return hash_->get_key_len() / 2; }

  class thread_type {
    inv_hash_storage_t* hash_;
    kMerUnitigInfo_type& info_;
    friend class jf_kmer_unitig_info;
    thread_type(inv_hash_storage_t* hash, kMerUnitigInfo_type& info) : hash_(hash), info_(info) { }

  public:
    void set(const mer_dna& mer, uint32_t number, uint32_t offset, uint32_t ori) {
      uint64_t id = 0, val;
      if(!hash_->get_val(mer[0], id, val, false))
        std::cerr << "Mer " << mer << " was not found in jf database\n";
      
      // There is no conflicts between threads because all k-mers are
      // unique in k-unitigs.
      auto unitig_info             = &info_[id];
      unitig_info->kUnitigNumber    = number;
      unitig_info->kUnitigOffset    = offset;
      unitig_info->kmerOriInKunitig = ori;
    }
  };
  thread_type thread() { return thread_type(hash_, kMerUnitigInfo_); }

  virtual const kMerUnitigInfoStruct* find(const mer_dna& mer) {
    uint64_t id, val;
    if(!hash_->get_val(mer[0], id, val, false))
      return 0;
    kMerUnitigInfoStruct* ret = &kMerUnitigInfo_[id];
    if(ret->kUnitigNumber == 0 && ret->kUnitigOffset == 0)
      return 0;
    return ret;
  }
};

class large_kmer_unitig_info : public mer_unitig_info_base {
  typedef multi_thread_skip_list_map<mer_dna, kMerUnitigInfoStruct> mer_info_map_type;
  mer_info_map_type kMerUnitigInfo_;
public:
  large_kmer_unitig_info() { }

  class thread_type {
    mer_info_map_type::thread thread_;
    friend class large_kmer_unitig_info;
    thread_type(mer_info_map_type& map) : thread_(map) { }

  public:
    void set(const mer_dna& mer, uint32_t number, uint32_t offset, uint32_t ori) {
      kMerUnitigInfoStruct info;
      info.kUnitigNumber    = number;
      info.kUnitigOffset    = offset;
      info.kmerOriInKunitig = ori;
      auto res = thread_.insert(std::make_pair(mer, info));
      if(!res.second)
        die << "kmer " << mer.to_str() << " already present in map";
    }
  };
  thread_type thread() { return thread_type(kMerUnitigInfo_); }

  virtual const kMerUnitigInfoStruct* find(const mer_dna& mer) {
    auto it = kMerUnitigInfo_.find(mer);
    if(it == kMerUnitigInfo_.end())
      return 0;
    return &it->second;
  }
};


ExpBuffer<uint64_t> kUnitigLengths;
int            readNumber; // TODO: Remove: it is not really used
int            longOutput;

template<typename mer_unitig_info_type>
class ReadKunitigs : public thread_exec {
  read_parser           unitig_parser;
  mer_unitig_info_type& mer_unitig_info;
  
public:
  ReadKunitigs(const char* kUnitigFile, mer_unitig_info_type& info, int nb_threads) :
    unitig_parser(kUnitigFile, nb_threads), mer_unitig_info(info) { }

  virtual void start(int th_id) {
    read_parser::stream unitig_stream(unitig_parser);
    mer_dna             fwd_mer;
    mer_dna             rev_mer;
    auto                thread_unitig_info = mer_unitig_info.thread();

    for( ; unitig_stream; ++unitig_stream) {
      // Parse header
      int kUnitigNumber, kUnitigLength;
      int fields_read = sscanf(unitig_stream->header, ">%d length:%d", &kUnitigNumber, &kUnitigLength);
      kUnitigLengths[kUnitigNumber] = kUnitigLength;
      if(fields_read != 2)
        die << "Fasta header of file does not match pattern '>UnitigiNumber length:UnitigLength'\n"
            << unitig_stream->header << "\n";
      if(unitig_stream->sequence.size() < (size_t)mer_dna::k())
        continue;
      char* cptr = unitig_stream->sequence;
      // Prime the k-mers
      for(int i = 0; i < (int)mer_dna::k() - 1; ++i, ++cptr) {
        fwd_mer.shift_left(*cptr);
        rev_mer.shift_right(mer_dna::complement(*cptr));
      }

      for( ; cptr < unitig_stream->sequence.end(); ++cptr) {
        fwd_mer.shift_left(*cptr);
        rev_mer.shift_right(mer_dna::complement(*cptr));
         
        const mer_dna* searchValue = &fwd_mer;
        unsigned char  ori         = 0;
        if(!(fwd_mer < rev_mer)) {
          searchValue = &rev_mer;
          ori         = 1;
        }
        thread_unitig_info.set(*searchValue, kUnitigNumber, cptr - unitig_stream->sequence.begin() - (mer_dna::k() - 1),
                               ori);
      }
    }
  }
};

/** Keep track whether or not a header has been written. On the first
    actual output (operator<<), the header is prepended. set_header()
    must be called to set/create the header. After we are done with a
    read, flush() is called, which adds a new line if necessary and
    endr (end of record) if asked to.
 */
class header_output {
  jflib::omstream& out_;
  charb            header_;
  bool             is_written_; // whether header has been written already
public:
  header_output(jflib::omstream& out) : out_(out), is_written_(false) { }

  void set_header(const char* name, size_t length) {
    sprintf(header_, "%s %ld", name, length);
    is_written_ = false;
  }
  void flush(bool end_of_record) {
    if(is_written_)
      out_ << "\n";
    if(end_of_record)
      out_ << jflib::endr;
  }

  template<typename T>
  jflib::omstream& operator<<(T& x) {
    if(!is_written_) {
      out_ << header_;
      is_written_ = true;
    }
    out_ << x;
    return out_;
  }
};

class ProcessReads : public thread_exec {
  read_parser           parser;
  mer_unitig_info_base* mer_unitig_info;
  jflib::o_multiplexer  multiplexer;

public:
  template<typename Iterator>
  ProcessReads(Iterator file_start, Iterator file_end, 
               mer_unitig_info_base* info, std::ostream& out_stream, int nb_threads) : 
    parser(file_start, file_end, nb_threads), 
    mer_unitig_info(info),
    multiplexer(&out_stream, 3 * nb_threads, 4096)
  {}

     virtual void start(int id) {
       read_parser::stream read_stream(parser);
       jflib::omstream     out(multiplexer);
       header_output       hout(out);
       char                read_prefix[3], prev_read_prefix[3];
       uint64_t            read_id = 0, prev_read_id = 0;
       memset(read_prefix, '\0', sizeof(read_prefix));
       
       for( ; read_stream; ++read_stream) {
         memcpy(prev_read_prefix, read_prefix, sizeof(read_prefix));
         prev_read_id = read_id;
         
         strtok(read_stream->header, " ");
	 if ((read_stream->header)[0] == '>')
	      sscanf(read_stream->header, ">%2s%ld", read_prefix, &read_id);
	 else
	      sscanf(read_stream->header, "@%2s%ld", read_prefix, &read_id);
         
         // Start new line & print read header if new read
         bool is_new_read = strcmp(read_prefix, prev_read_prefix) || read_id != prev_read_id;
         if(is_new_read) {
           // Keep mated read together in output.
           bool end_of_record =
             (read_id % 2 == 0) || (read_id != (prev_read_id + 1)) || strcmp(read_prefix, prev_read_prefix);
           hout.flush(end_of_record);
           hout.set_header(read_stream->header + 1, read_stream->sequence.size());
         }

         getMatchesForRead(read_stream->sequence, read_stream->sequence.end(),
                           read_stream->header + 1, mer_unitig_info, hout);
       }
       hout.flush(true);
     }
};

int readLastKUnitigNumber(const char* path) {
     int ret;
     FILE* infile = fopen (path, "r");
     if(!infile)
	  die << "Failed to open file '" << path << "'" << err::no;
     int fields_read = fscanf (infile, "%d", &ret);
     if(fields_read != 1)
	  die << "Failed to read the last k-unitig number from file '"
	      << path << "'" << err::no;
     fclose (infile);
     return ret;
}

int main(int argc, char *argv[])
{
     cmdline_parse  args(argc, argv);

     longOutput = args.long_flag;

     if(args.mer_arg == 0 || args.mer_arg > 10000)
       die << "Mer_len '" << args.mer_arg << "' is out of range";
     mer_dna::k(args.mer_arg);

     // Find out the last kUnitig number
     {
       int lastKUnitigNumber = readLastKUnitigNumber(args.numKUnitigsFile_arg);
       if(args.verbose_flag)
         std::cerr << "The largest kUnitigNumber is " <<  lastKUnitigNumber << "\n";
       kUnitigLengths.reserve(lastKUnitigNumber + 2);
     }

     if(args.verbose_flag)
       std::cerr << "mer length: " << mer_dna::k() << "\n";

     std::auto_ptr<mer_unitig_info_base> mer_unitig_info;
     if(args.jellyfishdb_given) {
       if(args.mer_arg > 31)
         die << "Mer size of " << args.mer_arg << " is incompatible with jellyfishdb switch (must have m <= 31).";
       auto jf_info = new jf_kmer_unitig_info(args.jellyfishdb_arg);
       mer_unitig_info.reset(jf_info);
       if(mer_dna::k() != jf_info->mer_len())
         die << "Mer length '" << mer_dna::k() << "' does not match jellyfish's mer length '" << jf_info->mer_len() << "'";
       ReadKunitigs<jf_kmer_unitig_info> KUnitigReader(args.kUnitigFile_arg, *jf_info, args.threads_arg);
       KUnitigReader.exec_join(args.threads_arg);
     } else {
       auto large_info = new large_kmer_unitig_info();
       mer_unitig_info.reset(large_info);
       ReadKunitigs<large_kmer_unitig_info> KUnitigReader(args.kUnitigFile_arg, *large_info, args.threads_arg);
       KUnitigReader.exec_join(args.threads_arg);
     }

     
     // Open output file and make sure it is deleted/closed on exit
     std::auto_ptr<std::ostream> out;
     if(args.gzip_flag)
       out.reset(new gzipstream(args.output_arg));
     else
       out.reset(new std::ofstream(args.output_arg));
     
     ProcessReads process_reads(args.readFiles_arg.begin(), args.readFiles_arg.end(),
                                mer_unitig_info.get(), *out, args.threads_arg);
     process_reads.exec_join(args.threads_arg);
     
     return (0);
}

void getMatchesForRead (const char *readBasesBegin, const char *readBasesEnd, 
                        const char *readName, mer_unitig_info_base* mer_unitig_info,
                        header_output& out)
{
     int readLength;
     mer_dna readKmer, readRevCompKmer;
     mer_dna readKmerTmp, readRevCompKmerTmp;
     unsigned char kUnitigOri, netOri, netOriHold, ori;
     int nextWithN;
     int jtmp, ktmp, tempIndex;
     int kUnitigNumber, kUnitigOffset;
     int isFirstRecord;
     int kUnitigOffsetOfFirstBaseInRead, kUnitigOffsetOfLastBaseInRead;
     int kUnitigOffsetOfFirstBaseInReadHold=0, kUnitigNumberHold=0;
     int intervalEnd=0;
     int ahg=0, bhg=0;
     int lastReadOffsetForNucmer=0;

#if 0 // only for debugging
     int intervalBegin=0;
     int firstKUnitigOffsetForNucmer=0, lastKUnitigOffsetForNucmer=0;
     int firstReadOffsetForNucmer=0;
#endif

     readLength = readBasesEnd - readBasesBegin;
     if((unsigned int)readLength < mer_dna::k())
       return;

     if (longOutput)
         out << "readNumber = " << readNumber 
             << " readLength = " << readLength << "\n";
       //fprintf (out, "readNumber = %d, readLength = %d\n", readNumber, readLength);
     //     readKmer = readRevCompKmer = 0;
     for (unsigned int j=0; j<mer_dna::k(); j++)
       readKmer.shift_left(readBasesBegin[j]);
     readRevCompKmer = readKmer;
     readRevCompKmer.reverse_complement();

     netOriHold = 2;
     nextWithN = readLength;
     for (int j=0; j<readLength; j++)
	  if (readBasesBegin[j] == 'N') {
	       nextWithN = j;
	       break; }
     for (int j=0, k=readLength-mer_dna::k(); k>=0; j++, k--) {
	  // Re-setting nextWithN if necessary
	  if (nextWithN < j)
	       for (jtmp=j; jtmp<readLength; jtmp++)
		    if (readBasesBegin[jtmp] == 'N') {
			 nextWithN = jtmp;
			 break; }
	  if (nextWithN < j)
	       nextWithN = readLength;
	  if (j>0) {
            readKmer.shift_left(readBasesBegin[j+mer_dna::k()-1]);
            readRevCompKmer.shift_right(mer_dna::complement(readBasesBegin[j+mer_dna::k()-1]));
          }

	  if (nextWithN-j < (int)mer_dna::k())
	       continue;
          const mer_dna* searchValue = &readKmer;
          ori = 0;
	  if (!(readKmer < readRevCompKmer)) {
	       searchValue = &readRevCompKmer;
	       ori = 1; }
          auto unitig_info = mer_unitig_info->find(*searchValue);
          if(!unitig_info)
            continue;
          kUnitigNumber = unitig_info->kUnitigNumber;
          kUnitigOffset = unitig_info->kUnitigOffset;
          kUnitigOri    = unitig_info->kmerOriInKunitig;

	  // now work with ori, j, and readNumber
	  isFirstRecord = 0;
	  if (kUnitigOri == ori)
	       netOri = 0;
	  else
	       netOri = 1;
	  // in the following kUnitigOffsetOfFirstBaseInRead
	  if (netOri == 0) // The read is 'F' in kUnitig
	       kUnitigOffsetOfFirstBaseInRead = kUnitigOffset - j;
	  else
	       kUnitigOffsetOfFirstBaseInRead = kUnitigOffset + j + mer_dna::k();
	  if ((kUnitigOffsetOfFirstBaseInReadHold != kUnitigOffsetOfFirstBaseInRead) ||
	      (kUnitigNumber != kUnitigNumberHold) ||
	      (netOri != netOriHold) ||
	      (j > intervalEnd)) {
	       if (netOriHold < 2) {
#if 0
		    if (longOutput)
			 fprintf (out, "Portion of read covered by the last k-unitig match ");
		    fprintf (out, "= (%d, %d)", intervalBegin, intervalEnd);
		    if (longOutput)
			 fprintf (out, "\nmyNucmerLine (ahg,bhg) =");
		    fprintf (out, " (%d,%d): %d %d", ahg, bhg, firstKUnitigOffsetForNucmer, lastKUnitigOffsetForNucmer);
		    if (netOriHold == 0)
			 fprintf (out, " %d %d", firstReadOffsetForNucmer, lastReadOffsetForNucmer);
		    else
			 fprintf (out, " %d %d", lastReadOffsetForNucmer, firstReadOffsetForNucmer);
		    fprintf (out, " %lu %d %d %s", kUnitigLengths[kUnitigNumberHold], readLength, kUnitigNumberHold, readName);
		    fprintf (out, "\n");
#else
                    out << " " << kUnitigNumberHold
                        << " " << ahg
                        << " " << ((netOriHold == 0) ? 'F' : 'R');
                    //		    fprintf (out, "%c %d %d %d %s\n", (netOriHold == 0) ? 'F' : 'R', ahg, readLength, kUnitigNumberHold, readName);
#endif
	       }
	       kUnitigOffsetOfFirstBaseInReadHold = kUnitigOffsetOfFirstBaseInRead;
	       kUnitigNumberHold = kUnitigNumber;
	       netOriHold = netOri;
#if 0
	       intervalBegin = j; 
#endif
	       isFirstRecord = 1; }
	  intervalEnd = j + mer_dna::k();
	  if (netOri == 0)
	       kUnitigOffsetOfLastBaseInRead = kUnitigOffsetOfFirstBaseInRead + readLength;
	  else
	       kUnitigOffsetOfLastBaseInRead = kUnitigOffsetOfFirstBaseInRead - readLength;
	  // Do Celera-style end calculations
	  // Use the k-unitig as read1 and the read as read2
	  if (netOri == 0) {
	       ahg = kUnitigOffsetOfFirstBaseInRead;
	       bhg = kUnitigOffsetOfLastBaseInRead - kUnitigLengths[kUnitigNumber]; }
	  else {
	       ahg = kUnitigOffsetOfLastBaseInRead;
	       bhg = kUnitigOffsetOfFirstBaseInRead - kUnitigLengths[kUnitigNumber]; }
	  ////////// BEGIN NUCMER END CALCULATIONS
	  if (ahg >= 0) {
#if 0
	       firstKUnitigOffsetForNucmer = ahg + 1;
#endif
	       if (netOri == 0) {
#if 0
		    firstReadOffsetForNucmer = 1;
#endif
               } else
		    lastReadOffsetForNucmer = readLength;
	  }
	  else {
#if 0
	       firstKUnitigOffsetForNucmer = 1;
#endif
	       if (netOri == 0) {
#if 0
		    firstReadOffsetForNucmer = 1 - ahg;
#endif
               } else
		    lastReadOffsetForNucmer = readLength + ahg;
	  }
	  if (bhg  >= 0) {
#if 0
            lastKUnitigOffsetForNucmer = kUnitigLengths[kUnitigNumber];
#endif
	       if (netOri == 0)
		    lastReadOffsetForNucmer = readLength - bhg;
#if 0
	       else
		    firstReadOffsetForNucmer = 1 + bhg;
#endif
	  }
	  else {
#if 0
	       lastKUnitigOffsetForNucmer = kUnitigLengths[kUnitigNumber] + bhg;
#endif
	       if (netOri == 0)
		    lastReadOffsetForNucmer = readLength;
#if 0
	       else
		    firstReadOffsetForNucmer = 1;
#endif
	  }
	  ////////// END NUCMER END CALCULATIONS
	  // Do we match at the end of the possible interval? (can jump)
	  if (isFirstRecord) {
	       jtmp = lastReadOffsetForNucmer - mer_dna::k();
	       ktmp = readLength - jtmp - mer_dna::k();
	       // If the desired k-mer has an 'N', don't bother
	       if (nextWithN < jtmp) {
                   for (tempIndex=jtmp; tempIndex<jtmp+(int)mer_dna::k(); tempIndex++)
			 if (readBasesBegin[tempIndex] == 'N')
			      break; }
	       else
		    tempIndex = nextWithN;
	       if (tempIndex-jtmp < (int)mer_dna::k()) // It had an 'N'
		    continue;
	       readKmerTmp = readKmer;
	       readRevCompKmerTmp = readRevCompKmer;
	       tempIndex = j + mer_dna::k() - jtmp;
	       if (tempIndex < 0)
		    tempIndex = 0;
	       for (tempIndex=0; tempIndex<(int)mer_dna::k(); tempIndex++) {
                 readKmerTmp.shift_left(readBasesBegin[jtmp+tempIndex]);
                 readRevCompKmerTmp.shift_right(mer_dna::complement(readBasesBegin[jtmp+tempIndex]));
               }
	 
               searchValue = &readKmerTmp;
               ori = 0;
	       if (!(readKmerTmp < readRevCompKmerTmp)) {
		    searchValue = &readRevCompKmerTmp;
		    ori = 1; }
               auto unitig_info = mer_unitig_info->find(*searchValue);
               if(!unitig_info)
                 continue;
               kUnitigNumber = unitig_info->kUnitigNumber;
               kUnitigOffset = unitig_info->kUnitigOffset;
               kUnitigOri    = unitig_info->kmerOriInKunitig;

	       // now work with ori, j, and readNumber
	       if (kUnitigOri == ori)
		    netOri = 0;
	       else
		    netOri = 1;
	       // in the following kUnitigOffsetOfFirstBaseInRead
	       if (netOri == 0) // The read is 'F' in kUnitig
		    kUnitigOffsetOfFirstBaseInRead = kUnitigOffset - jtmp;
	       else
		    kUnitigOffsetOfFirstBaseInRead = kUnitigOffset + jtmp + mer_dna::k();
	       if ((kUnitigOffsetOfFirstBaseInReadHold == kUnitigOffsetOfFirstBaseInRead) &&
		   (kUnitigNumber == kUnitigNumberHold) &&
		   (netOri == netOriHold)) {
		    j = jtmp;
		    k = ktmp; 
		    readKmer = readKmerTmp;
		    readRevCompKmer = readRevCompKmerTmp;
		    intervalEnd = j + mer_dna::k(); }
	  }
     }
     //	    break; // For debugging
     if (netOriHold < 2) {
#if 0
	  if (longOutput)
	       fprintf (out, "Portion of read covered by the last k-unitig match ");
	  fprintf (out, "= (%d, %d)", intervalBegin, intervalEnd);
	  if (longOutput)
	       fprintf (out, "\nmyNucmerLine (ahg,bhg) =");
	  fprintf (out, " (%d,%d): %d %d", ahg, bhg, firstKUnitigOffsetForNucmer, lastKUnitigOffsetForNucmer);
	  if (netOriHold == 0)
	       fprintf (out, " %d %d", firstReadOffsetForNucmer, lastReadOffsetForNucmer);
	  else
	       fprintf (out, " %d %d", lastReadOffsetForNucmer, firstReadOffsetForNucmer);
	  fprintf (out, " %lu %d %d %s", kUnitigLengths[kUnitigNumberHold], readLength, kUnitigNumberHold, readName);
	  fprintf (out, "\n");
#else
          //	  fprintf (out, "%c %d %d %d %s\n", (netOriHold == 0) ? 'F' : 'R', ahg, readLength, kUnitigNumberHold, readName);
          out << " " << kUnitigNumberHold
              << " " << ahg
              << " " << ((netOriHold == 0) ? 'F' : 'R');
#endif
     }
}

