#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include <unittests/test_main.hpp>

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  test_main args(argc, argv);
  
  unsigned int seed;
  if(args.seed_given) {
    seed = args.seed_arg;
  } else {
    std::ifstream urandom("/dev/urandom");
    urandom.read((char*)&seed, sizeof(seed));
    if(!urandom.good()) {
      std::cerr << "Failed to read random seed" << std::endl;
      return 1;
    }
  }

  std::cout << "Using random seed " << seed << std::endl;
  srandom(seed);
  
  return RUN_ALL_TESTS();
}
