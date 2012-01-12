#ifndef __SKIP_LIST_SET_HPP
#define __SKIP_LIST_SET_HPP

/* An ordered set class based on the skip list data structure, instead
   of the more usual balanced binary tree. The performance is similar
   to a balanced tree, i.e. log(n) operations for insertion and
   search, but in average instead of worst case.  

   The interface is almost a drop-in replacement for std::set, with
   the exception that the iterators forward iterators instead of
   bidirectional iterators. And consequently there is no reverse
   operators (no rbegin() and rend()).

   As it stands, there is no allocator parameter.
 */

#include <assert.h>
#include <iterator>
#include <algorithm>

// Count trailing zeros
template<typename T>
int ctz(T x);
template<>
int ctz<unsigned int>(unsigned int x) { return __builtin_ctz(x); }
template<>
int ctz<unsigned long>(unsigned long x) { return __builtin_ctzl(x); }
template<>
int ctz<unsigned long long>(unsigned long long x) { return __builtin_ctzll(x); }

// Random function
struct std_random {
  typedef unsigned long rand_type;
  rand_type operator()() const { return (rand_type)random(); }
};

// Return the height of a tower of pointer. Specialized for 2 and 4.
template<typename Random, int p>
struct random_height;
template<typename Random>
struct random_height<Random, 2> {
  Random rng;
  int operator()() const {
    typename Random::rand_type x = rng();
    return (x == 0 ? 8*sizeof(typename Random::rand_type) : ctz(x)) + 1;
  }
  random_height(const Random& rng_ = Random()) : rng(rng_) { }
};
template<typename Random>
struct random_height<Random, 4> {
  Random rng;
  int operator()() const {
    typename Random::rand_type x = rng();
    return (x == 0 ? 4 * sizeof(typename Random::rand_type) : ctz(x) >> 1) + 1;
  }
  random_height(const Random& rng_ = Random()) : rng(rng_) { }
};

// Set based on a skip list
template <typename Key, typename Compare = std::less<Key>, int p_ = 4, typename Random = std_random>
class skip_list_set {
  struct node {
    Key   k;
    int   height;
    node* tower[];
  };
  node**                    heads_;
  int                       max_height_;
  int                       cur_height_; // 0 < cur_height <= max_height_
  size_t                    size_;
  Compare                   comp_;
  random_height<Random, p_> rh_;

public:
  typedef Key                                   key_type;
  typedef Key                                   value_type;
  typedef Compare                               key_compare;
  typedef Compare                               value_compare;
  typedef Key&                                  reference;
  typedef const Key&                            const_reference;
  typedef size_t                                size_type;
  typedef ssize_t                               difference_type;
  typedef Key*                                  pointer;
  typedef const Key*                            const_pointer;
  static const int p = p_;

  // This template stuff for defining the iterator and const_iterator
  // is barely worth it. It is hardly shorter and more confusing.
  template<typename ref_type, typename ptr_type>
  class node_iterator : public std::iterator<std::forward_iterator_tag, key_type> {
    node* item;
    friend class skip_list_set;
    node_iterator(node* item_) : item(item_) { }
  public:
    node_iterator() : item(0) { }
    node_iterator(const node_iterator& rhs) : item(rhs.item) { }

    node_iterator& operator=(const node_iterator& rhs) {
      item = rhs.item;
      return *this;
    }    
    bool operator==(const node_iterator& rhs) const { return item == rhs.item; }
    bool operator!=(const node_iterator& rhs) const { return item != rhs.item; }
    ref_type operator*() { return item->k; }
    ptr_type operator->() { return &item->k; }
    node_iterator& operator++() {
      item = item->tower[0];
      return *this;
    }
    node_iterator operator++(int) {
      node_iterator c(*this);
      item = item->tower[0];
      return c;
    }
  };
  typedef node_iterator<reference, pointer> iterator;
  class const_iterator : public node_iterator<const_reference, const_pointer> {
    typedef node_iterator<const_reference, const_pointer> super;
    friend class skip_list_set;
    const_iterator(node* item_) : super(item_) { }
  public:
    const_iterator() : super(0) { }
    const_iterator(const const_iterator& rhs) : super(rhs.item) { }
    const_iterator(const iterator& rhs) : super(rhs.item) { }
    using super::operator==;
    using super::operator!=;
    bool operator==(const iterator& rhs) const { return super::item == rhs.item; }
    bool operator!=(const iterator& rhs) const { return super::item != rhs.item; }
  };

  explicit skip_list_set(const Compare& comp = Compare(), 
                         const Random& rand = Random()) :
    heads_(new node*[10]), max_height_(10), cur_height_(1), size_(0),
    comp_(comp), rh_(rand)
  {
    for(int i = 0; i < max_height_; ++i)
      heads_[i] = 0;
    
  }
                         
  explicit skip_list_set(int max_height,
                         const Compare& comp = Compare(), 
                         const Random& rand = Random()) :
    heads_(new node*[max_height]),
    max_height_(max_height), cur_height_(1), size_(0),
    comp_(comp), rh_(rand)
  {
    for(int i = 0; i < max_height_; ++i)
      heads_[i] = 0;
  }
  template<class InputIterator>
  skip_list_set(int max_height, InputIterator first, InputIterator last,
                const Compare& comp = Compare(), 
                const Random& rand = Random()) :
    heads_(new node*[max_height]),
    max_height_(max_height), cur_height_(1), size_(0),
    comp_(comp), rh_(rand)
  {
    for(int i = 0; i < max_height_; ++i)
      heads_[i] = 0;
    insert(first, last);
  }    
  skip_list_set(const skip_list_set& rhs) :
    heads_(new node*[rhs.max_height_]),
    max_height_(rhs.max_height_), cur_height_(1), size_(0)
  {
    for(int i = 0; i < max_height_; ++i)
      heads_[i] = 0;
    insert(rhs.begin(), rhs.end());
  }
  virtual ~skip_list_set() {
    clear();
    delete [] heads_;
  }

  skip_list_set& operator=(skip_list_set rhs) {
    swap(rhs);
    return *this;
  }

  void swap(skip_list_set& rhs) {
    std::swap(heads_, rhs.heads_);
    std::swap(max_height_, rhs.max_height_);
    std::swap(cur_height_, rhs.cur_height_);
    std::swap(size_, rhs.size_);
    std::swap(comp_, rhs.comp_);
    std::swap(rh_, rhs.rh_);
  }

  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  iterator begin() { return iterator(heads_[0]); }
  const_iterator begin() const { return const_iterator(heads_[0]); }
  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }

  void clear() {
    node* cnode = heads_[0];
    while(cnode) {
      node* nnode = cnode->tower[0];
      delete cnode;
      cnode = nnode;
    }
    memset(heads_, '\0', sizeof(node*) * max_height_);
    size_ = 0;
  }

  size_type count(const key_type& x) const {
    node **path[max_height_];
    if(find_node_path(x, path))
      return 1;
    return 0;
  }

  std::pair<iterator, iterator> equal_range(const key_type& x) const {
    node **path[max_height_];
    node *n = find_node_path(x, path);
    if(n)
      return std::make_pair(iterator(n), ++iterator(n));
    return std::make_pair(iterator(*path[0]), iterator(*path[0]));
  }

  size_type erase(const key_type& x) {
    node **path[max_height_];
    node *n = find_node_path(x, path);
    if(!n)
      return 0;
    for(int i = 0; i < n->height; ++i)
      *path[i] = n->tower[i];
    delete n;
    --size_;
    return 1;
  }

  iterator find(const key_type& x) const {
    node** path[max_height_];
    return iterator(find_node_path(x, path));
  }

  key_compare key_comp() const { return comp_; }
  value_compare value_comp() const { return comp_; }

  iterator lower_bound(const key_type& x) const {
    node** path[max_height_];
    find_node_path(x, path);
    return iterator(*path[0]);
  }

  iterator upper_bound(const key_type& x) const {
    node** path[max_height_];
    node* n =find_node_path(x, path);
    if(n)
      return ++iterator(n);
    return iterator(*path[0]);
  }

  size_type max_size() const {
    size_type res = 1;
    size_type pp  = p;
    int       n   = max_height_;

    while(n) {
      if(n & 0x1)
        res *= pp;
      pp *= pp;
      n >>= 1;
    }
    return res;
  }

  std::pair<iterator, bool> insert(const value_type& x) {
    node** path[max_height_];
    node*  n = find_node_path(x, path);
    if(n)
      return std::make_pair(iterator(n), false);
    n = new_node(x);
    for(int i = 0; i < std::min(n->height, cur_height_); ++i) {
      n->tower[i] = *path[i];
      *path[i] = n;
    }
    if(n->height > cur_height_) {
      n->tower[cur_height_] = 0;
      heads_[cur_height_++] = n;
    }
    ++size_;
    return std::make_pair(iterator(n), true);
  }

  // Unlike std::set, this provides no advantage at this point.
  iterator insert(iterator position, const value_type& x) {
    return insert(x).first;
  }

  template<class InputIterator>
  void insert(InputIterator first, InputIterator last) {
    for( ; first != last; ++first)
      insert(*first);
  }
  
private:
  // Find the node path. I.e. the addresses of the pointers that leads
  // to the largest element less than x. If a node with a value
  // equivalent to x is found, a pointer to the node is returned and
  // path has an undefined value. Otherwise, 0 is returned and the
  // path array is initialized.
  node* find_node_path(const value_type& x, node*** path) const {
    node*  cnode = 0;
    for(int i = cur_height_ - 1; i >= 0; --i) {
      node** cptr = cnode ? &cnode->tower[i] : &heads_[i];
      while(*cptr) {
        if(!comp_((*cptr)->k, x))
          break;
        cnode = *cptr;
        cptr  = &(cnode->tower[i]);
      }
      path[i] = cptr;
    }
    // Check if we found a node equal to x. If so, return it.
    return *path[0] && !comp_(x,(*path[0])->k) ? *path[0] : 0 ;
  }

  // Allocate a new node. Does raw memory allocation of a node with
  // enough space for the tower. Then in place copy construction of
  // the key from x.
  node* new_node(const value_type& x) const {
    int height  = std::min(max_height_, std::min(cur_height_ + 1, rh_()));
    node* res   = (node*)operator new(sizeof(node) + height * sizeof(node*));
    res->height = height;
    new ((void*)&res->k) value_type(x);
    return res;
  }

  // // Debugging routines.
  // void print_lists(std::ostream &os = std::cout) const {
  //   for(int i = max_height_ - 1; i >= 0; --i) {
  //     node* cnode = heads_[i];
  //     os << i << ": " << (void*)&heads_[i];
  //     do {
  //       os << " " << (void*)cnode;
  //       if(cnode) {
  //         os << ":" << cnode->k;
  //         cnode = cnode->tower[i];
  //       }
  //     } while(cnode);
  //     os << "\n";
  //   }
  //   os << std::flush;
  // }

  // void check_path(const value_type& x, node*** path) const {
  //   for(int i = cur_height_ - 1; i >= 0; --i) {
  //     std::cout << "check " << i << ": " << (void*)path[i] << ":"
  //               << (path[i] ? (void*)*path[i] : (void*)0) << std::endl;
  //   }
  // }
};


#endif /* __SKIP_LIST_SET_HPP */
