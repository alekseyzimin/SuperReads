#include<heap.hpp>
#include<gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>

TEST(Heap, MaxHeap) {
  heap<int>::max h;
  std::vector<int> elts;

  for(int i = 0; i < 10; ++i)
    elts.push_back(i);
  std::random_shuffle(elts.begin(), elts.end());

  int max = std::numeric_limits<int>::min();
  for(int i = 0; i < (int)elts.size(); ++i) {
    max = std::max(max, elts[i]);
    h.push(elts[i]);
    EXPECT_EQ(i+1, (int)h.size());
    EXPECT_EQ(max, h.peek());
  }
  for(int i = 0; i < (int)elts.size(); ++i)
    EXPECT_EQ((int)elts.size() - i - 1, h.pop());
}

TEST(Heap, MinHeap) {
  heap<int>::min h;
  std::vector<int> elts;

  for(int i = 0; i < 10; ++i)
    elts.push_back(i);
  std::random_shuffle(elts.begin(), elts.end());

  int min = std::numeric_limits<int>::max();
  for(int i = 0; i < (int)elts.size(); ++i) {
    min = std::min(min, elts[i]);
    h.push(elts[i]);
    EXPECT_EQ(i+1, (int)h.size());
    EXPECT_EQ(min, h.peek());
  }

  for(int i = 0; i < (int)elts.size(); ++i)
    EXPECT_EQ(i, h.pop());
}

TEST(Heap, Heapify) {
  std::vector<int> elts;

  for(int i = 0; i < 10; ++i)
    elts.push_back(i);
  std::random_shuffle(elts.begin(), elts.end());

  heap<int>::min h(elts.begin(), elts.end());
  EXPECT_EQ(elts.size(), h.size());
  for(int i = 0; i < (int)elts.size(); ++i)
    EXPECT_EQ(i, h.pop());
}

TEST(Heap, Copy) {
  std::vector<int> elts;

  for(int i = 0; i < 10; ++i)
    elts.push_back(i);
  std::random_shuffle(elts.begin(), elts.end());

  heap<int>::min h1(elts.begin(), elts.end());
  heap<int>::min h2;

  EXPECT_EQ(elts.size(), h1.size());
  EXPECT_EQ((size_t)0, h2.size());

  std::swap(h1, h2);
  EXPECT_EQ(elts.size(), h2.size());
  EXPECT_EQ((size_t)0, h1.size());

  heap<int>::min h3(h2);
  EXPECT_EQ(elts.size(), h2.size());
  EXPECT_EQ(elts.size(), h3.size());
  
  heap<int>::min h4;
  h4 = h2;
  EXPECT_EQ(elts.size(), h2.size());
  EXPECT_EQ(elts.size(), h4.size());
}
