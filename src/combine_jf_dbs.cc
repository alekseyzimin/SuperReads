#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/compacted_hash.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/misc.hpp>
#include <src/combine_jf_dbs.hpp>

int main(int argc, char *argv[])
{
  combine_jf_dbs     args(argc, argv);
  SquareBinaryMatrix hash_matrix;
  unsigned int       reprobe;
  unsigned int       klen;
  unsigned int       N = args.db_jf_arg.size();
  inv_hash_storage_t *ary;

  if(args.verbose_flag)
    std::cerr << "Merging " << N << " databases\n";

  // Get information from the first database. Other database must
  // match.
  auto db_it = args.db_jf_arg.begin();
  {
    mapped_file dbf(*db_it);
    dbf.will_unmap(true);
    dbf.sequential().will_need();
    char type[8];
    memcpy(type, dbf.base(), sizeof(type));
    if(!strncmp(type, jellyfish::raw_hash::file_type, sizeof(type))) {
      raw_inv_hash_query_t hash(dbf);
      klen        = hash.get_key_len();
      reprobe     = hash.get_max_reprobe();
      hash_matrix = hash.get_hash_matrix();
      if(args.verbose_flag)
        std::cerr << "Key len: " << klen << " reprobe:" << reprobe << "\n";
      ary = new inv_hash_storage_t(hash.get_size(), klen, 
                                   hash.get_val_len() + ceilLog2(N-1),
                                   reprobe, jellyfish::quadratic_reprobes);
      ary->set_matrix(hash_matrix);
      if(args.verbose_flag)
        std::cerr << "Loading database: " << *db_it << "\n";
      auto it = hash.get_iterator();
      while(it.next())
        if(!ary->map(it.get_key(), it.get_val() * N))
          die << "Output hash is full";
    } else {
      die << "Invalid file type '" << err::substr(type, sizeof(type)) << "'.";
    }
  }

  ++db_it;
  for(uint64_t nb = 1; db_it != args.db_jf_arg.end(); ++db_it, ++nb) {
    mapped_file dbf(*db_it);
    dbf.will_unmap(true);
    dbf.sequential().will_need();
    char type[8];
    memcpy(type, dbf.base(), sizeof(type));
    if(!strncmp(type, jellyfish::raw_hash::file_type, sizeof(type))) {
      raw_inv_hash_query_t hash(dbf);
      if(klen != hash.get_key_len())
        die << "Invalid key len for '" << *db_it << "': " << hash.get_key_len()
            << " expected " << klen;
      if(reprobe != hash.get_max_reprobe())
        die << "Invalid max reprobe for '" << *db_it << "': " << hash.get_max_reprobe()
            << " expected " << reprobe;
      if(hash_matrix != hash.get_hash_matrix())
        die << "Invalid hash matrix for '" << *db_it << "'";
      if(args.verbose_flag)
        std::cerr << "Loading database: " << *db_it << "\n";
      auto it = hash.get_iterator();
      while(it.next())
        if(!ary->map(it.get_key(), it.get_val() * N + nb))
          die << "Output hash is full";
    } else {
      die << "Invalid file type '" << err::substr(type, sizeof(type)) << "'.";
    }    
  }
  
  if(args.verbose_flag)
    std::cerr << "Writing result to '" << args.output_arg << "'\n";
  raw_inv_hash_dumper_t dumper((uint_t)4, args.output_arg, 10000000, ary);
  dumper.dump();

  return 0;
}
