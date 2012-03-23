#ifndef __MULTI_THREAD_SKIP_LIST_MAP_HPP__
#define __MULTI_THREAD_SKIP_LIST_MAP_HPP__

#include <multi_thread_skip_list_set.hpp>

/** Multi-threaded lock free map based on skip list. See
    multi_thread_skip_list_set.hpp for more details.
*/

template<typename Key, typename T, class Compare = std::less<Key>,
         int p_ = 4, typename Random = xor_random>
class multi_thread_skip_list_map :
  public multi_thread_skip_list_set<std::pair<const Key, T>, first_comp<std::pair<const Key,T>, Compare>, p_, Random>
{
  typedef multi_thread_skip_list_set<std::pair<const Key, T>, first_comp<std::pair<const Key,T>, Compare>, p_, Random> super;
public:
  typedef Key                             key_type;
  typedef T                               mapped_type;
  typedef std::pair<const                 Key, T> value_type;
  typedef Compare                         key_compare;
  typedef first_comp<value_type, Compare> value_compare;
  typedef value_type&                     reference;
  typedef const value_type&               const_reference;
  typedef typename super::iterator        iterator;
  typedef typename super::const_iterator  const_iterator;
  typedef typename super::size_type       size_type;
  typedef typename super::difference_type difference_type;
  typedef typename super::pointer         pointer;
  typedef typename super::const_pointer   const_pointer;

  explicit multi_thread_skip_list_map(int max_height = 10,
                                      const Compare& comp = Compare(),
                                      const Random& rand = Random()) :
    super(max_height, value_compare(comp), rand)
  { }

  virtual ~multi_thread_skip_list_map() { }
  
  class thread : public super::thread {
  public:
    thread(multi_thread_skip_list_map& map) : super::thread(map) { }
    virtual ~thread() { }

    mapped_type& operator[](const key_type& x) {
      return (*((this->insert(std::make_pair(x,T()))).first)).second;
    }
  };
};

#endif /* __MULTI_THREAD_SKIP_LIST_MAP_HPP__ */

