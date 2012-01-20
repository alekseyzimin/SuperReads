#ifndef __HEAP_HPP__
#define __HEAP_HPP__

#include <vector>
#include <algorithm>
#include <functional>

template<typename T, typename Pr = std::less<T>, typename C = std::vector<T> >
class basic_heap {
public:
  typedef T                          value_type;
  typedef typename C::iterator       iterator;
  typedef typename C::const_iterator const_iterator;

  basic_heap() { }
  template<typename ForwardIterator>
  basic_heap(ForwardIterator begin, ForwardIterator end) : ary(begin, end) {
    make_heap(ary.begin(), ary.end(), comp);
  }
  virtual ~basic_heap() { }

  void swap(basic_heap &rhs) {
    std::swap(ary, rhs.ary);
    std::swap(comp, rhs.comp);
  }
  basic_heap & operator=(basic_heap rhs) {
    swap(rhs);
    return *this;
  }
  
  void push(const T &e) { 
    ary.push_back(e);
    push_heap(ary.begin(), ary.end(), comp);
  }
  T &pop() {
    pop_heap(ary.begin(), ary.end(), comp);
    T &tmp = ary.back();
    ary.pop_back();
    return tmp;
  }
  const T &peek() const { return ary.front(); }
  size_t size() const { return ary.size(); }
  bool empty() const { return ary.empty(); }

  iterator begin() { return ary.begin(); }
  const_iterator begin() const { return ary.begin(); }
  iterator end() { return ary.end(); }
  const_iterator end() const { return ary.end(); }
  void clear() { ary.clear(); }

protected:
  C  ary;
  Pr comp;
};

namespace std {
  template<typename T, typename Pr, typename C>
  void swap(basic_heap<T, Pr, C> &h1, basic_heap<T, Pr, C> &h2) {
    h1.swap(h2);
  } 
}

/* This is what we would like to do: do a typedef template. Alas,
   C++03 does not support it. The following two lines will not
   compile. */
// template<typename T, typename C = std::vector<T> >
// typedef basic_heap<T, std::less<T>, C> min_heap;

/* To work around it, define a struct template. Typedefs in this
   struct can depend on the template parameter of the struct. This
   defines two types: heap::min and heap::max. Close to what we
   originally want. */
template<typename T, typename C = std::vector<T> >
struct heap {
  typedef basic_heap<T, std::less<T>, C> max;
  typedef basic_heap<T, std::greater<T>, C> min;
};

#endif /* __HEAP_HPP__ */
