#include <gtest/gtest.h>
#include <unittests/misc.hpp>
#include <multi_thread_skip_list_map.hpp>

#include <string>
#include <vector>
#include <map>

namespace {
  typedef multi_thread_skip_list_map<int, int> map_type;

  struct thread_insert_data {
    map_type                 map; // The set to add to
    jflib::atomic_field<int> ids; // Counter to get a thread id
    std::vector<int>         v;   // The data to insert

    static const int per_th = 10000; // Nb elements per thread
    thread_insert_data(int nb_threads) :
      map(10, std::less<int>(), xor_random(random())), ids(-1)
    { 
      const int n = per_th * nb_threads;
      for(int i = 0; i < n; ++i)
        v.push_back(random() % n);
    }
  };

  void* insert_from_thread(void* d) {
    auto data = (thread_insert_data*)d;
    map_type::thread th(data->map);

    int tid = (data->ids += 1);
    for(int i = tid * data->per_th; i < (tid+1) * data->per_th; ++i)
      th[data->v[i]] += data->v[i];

    return 0;
  }

  TEST(MTSkipListMap, InsertManyThreads) {
    static const int   nb_threads = 5;
    thread_insert_data data(nb_threads);

    pdo(nb_threads, insert_from_thread, (void*)&data);
    
    // Do the same single threads into a set and check the statistics
    std::map<int, int> std_map;
    for(auto it = data.v.begin(); it != data.v.end(); ++it)
      std_map[*it] += *it;

    EXPECT_EQ(std_map.size(), data.map.size());

    auto it = data.map.begin();
    for(auto std_it = std_map.begin(); std_it != std_map.end(); ++std_it, ++it) {
      EXPECT_NE(data.map.end(), it);
      EXPECT_EQ(std_it->first, it->first);
      EXPECT_EQ(std_it->second, it->second);
    }
    EXPECT_EQ(data.map.end(), it);
  }

} // namespace
