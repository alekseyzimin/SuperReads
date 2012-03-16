#ifndef __EXP_VECTOR_HPP__
#define __EXP_VECTOR_HPP__

#include <memory>
#include <vector>

template <class T, class Allocator = std::allocator<T> >
class exp_vector : public std::vector<T, Allocator> {
  typedef std::vector<T, Allocator> super;
public:
  typedef typename super::reference              reference;
  typedef typename super::const_reference        const_reference;
  typedef typename super::iterator               iterator;
  typedef typename super::const_iterator         const_iterator;
  typedef typename super::size_type              size_type;
  typedef typename super::difference_type        difference_type;
  typedef typename super::value_type             value_type;
  typedef typename super::allocator_type         allocator_type;
  typedef typename super::pointer                pointer;
  typedef typename super::const_pointer          const_pointer;
  typedef typename super::reverse_iterator       reverse_iterator;
  typedef typename super::const_reverse_iterator const_reverse_iterator;

  explicit exp_vector(const Allocator& a = Allocator()) : super(a) { }
  explicit exp_vector (size_type n, const T& value = T(), const Allocator& a = Allocator()) :
    super(n, value, a) { }
  template <class InputIterator>
  exp_vector ( InputIterator first, InputIterator last, const Allocator& a = Allocator() ) :
    super(first, last, a) { }
  exp_vector ( const std::vector<T,Allocator>& x ) :
    super(x) { }

  virtual ~exp_vector() { }

  reference operator[] ( size_type n ) {
    if(n >= super::size())
      super::resize(n + 1);
    return super::operator[](n);
  }
  const_reference operator[] ( size_type n ) const {
    if(n >= super::size())
      super::resize(n + 1);
    return super::operator[](n);
  }
};

#endif /* __EXP_VECTOR_HPP__ */
