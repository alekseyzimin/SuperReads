#include <unittests/misc.hpp>

void pdo(unsigned int n, void* (*f)(void*), void *data) {
  pthread_t tids[n];
  for(unsigned int i = 0; i < n; ++i) {
    int e = pthread_create(&tids[i], 0, f, data);
    if(e)
      throw std::runtime_error(strerror(e));
  }
  for(unsigned int i = 0; i < n; ++i) {
    int e = pthread_join(tids[i], 0);
    if(e)
      throw std::runtime_error(strerror(e));
  }
}
