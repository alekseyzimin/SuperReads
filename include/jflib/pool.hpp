#ifndef _JFLIB_POOL_H_
#define _JFLIB_POOL_H_

#include <jflib/circular_buffer.hpp>
#include <jflib/locks_pthread.hpp>
#include <vector>
#include <stdexcept>

namespace jflib {
  /** Fixed size pool of objects with lock-free access. Conceptually,
   * there are two sides: A and B. There are 2 FIFOs in opposite
   * direction between A and B. One "get" an element from one side
   * (say A) and then "release" it. This correspond to dequeue from
   * the queue B->A then enqueueing into A->B. To get an element can
   * block if the queue is empty but to release never does.
   *
   * The class itself manages 'size' elements in the queues and are
   * initialized with 'value'. Initially, all the elements are in the
   * B->A queue.
   *
   * The elements of type T in the pool should not be accessed
   * unless one has the obtained the elt object and has not released
   * (or destructed) it. 
   *
   * WARNING: The iterator on the elements does not check any of this
   * and is not thread-safe.
   */
  template<typename T, typename CV = locks::cond>
  class pool {
    class                           side;

  public:
    typedef typename std::vector<T> Tvec;
    
    pool(size_t size, const T &value = T()) 
      : elts(size, value), B2A(size, &elts), A2B(size, &elts)
    { 
      B2A._other = &A2B;
      A2B._other = &B2A;
      for(uint32_t i = 0; i != elts.size(); ++i)
        B2A._fifo.enqueue(i);
    }
    virtual ~pool() { }

    /** Get the side A. Used to get an element from side A. */
    side &get_A() { return B2A; }
    /** Get the side B. Used to get an element from side B. */
    side &get_B() { return A2B; }
    void close_A_to_B() { A2B._fifo.close(); A2B.signal(true); }
    void close_B_to_A() { B2A._fifo.close(); B2A.signal(true); }
    bool is_closed_A_to_B() const { return A2B._fifo.is_closed(); }
    bool is_closed_B_to_A() const { return B2A._fifo.is_closed(); }

    /** Iterators on the elements. Unlike other 
     */
    typename Tvec::iterator begin() { return elts.begin(); }
    typename Tvec::iterator end() { return elts.end(); }

    /** A wrapper around an element of type T. The element can be
     * obtained with operator* or operator->. release() is called by
     * the destructor to requeue the element in to the double
     * fifo. When the double fifo is empty and closed, the element
     * obtained is_empty() method returns true.
     */
    class elt {
    public:
      elt() : _i(cbT::guard), _v(0), _s(0) { }
      elt(side &s) : _i(s.get()), _v(s[_i]), _s(s._other) { }
      ~elt() { release(); }
      elt &operator=(side &s) {
        release();
        _i = s.get();
        _v = s[_i];
        _s = s._other;
        return *this;
      }

      void release() { 
        if(_v)
          _s->release(_i);
        _v = 0;
      }
      bool is_empty() { return _v == 0; }
      T &operator*() { return *_v; }
      T *operator->() { return _v; }

      friend class pool;
    private:
      elt(const elt &rhs) { }
      elt &operator=(const elt &rhs) { }

      uint32_t  _i;             // Index of stored value
      T        *_v;             // Stored value
      side     *_s;             // Side to release to
    };
    static const elt closed;

  private:
    typedef circular_buffer<uint32_t> cbT;
    /** A circular buffer with a conditional variable to wait in the
     * empty event and a pointer to the other direction circular
     * buffer. Every method is private and are accessible only by its
     * friends :)
     */
    class side {
    private:
      friend class pool;
      friend class elt;
      enum State { NONE, WAITING, CLOSED };
      side(size_t size, Tvec *elts) : 
        _fifo(2*size), _state(NONE), _other(0), _elts(elts) { }

      uint32_t get();
      T *operator[](uint32_t i);
      void release(uint32_t i);
      void signal(bool force = false);

      cbT   _fifo;
      CV    _cond;
      State _state;
      side *_other;
      Tvec *_elts;
    };

    Tvec elts;
    side B2A;                     // Manages queue from B->A
    side A2B;                     // Manages queue from A->B
  };
}

template<typename T, typename CV>
uint32_t jflib::pool<T, CV>::side::get() {
  bool     last_attempt = false;
  uint32_t res          = _fifo.dequeue();
  while(res == cbT::guard) {
    _cond.lock();

    switch(a_get(_state)) {
    case CLOSED:
      if(last_attempt) {
        _cond.unlock();
        return cbT::guard;
      } else {
        last_attempt = true;
        break;
      }
    case NONE:
      a_set(_state, WAITING);
      break;
    case WAITING:
      break;
    }
    res = _fifo.dequeue();
    if(res == cbT::guard) {
      if(last_attempt) {
        _cond.unlock();
        break;
      }
    } else {
      _cond.unlock();
      break;
    }
    do {
      _cond.timedwait(5);
    } while(a_get(_state) == WAITING);
    _cond.unlock();
  }

  return res;
}

template<typename T, typename CV>
T * jflib::pool<T, CV>::side::operator[](uint32_t i) {
  if(i == cbT::guard)
    return 0;
  return &(*_elts)[i];
}

template<typename T, typename CV>
void jflib::pool<T, CV>::side::release(uint32_t i) {
  while(!_fifo.enqueue(i)) ;
  signal();
}

template<typename T, typename CV>
void jflib::pool<T, CV>::side::signal(bool close) {
  if(a_get(_state) != NONE || close) {
    _cond.lock();
    a_set(_state, close ? CLOSED : NONE);
    _cond.broadcast();
    _cond.unlock();
  }
}

#endif /* _JFLIB_POOL_H_ */
