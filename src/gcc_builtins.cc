#include <gcc_builtins.hpp>

template<>
int ctz<unsigned int>(unsigned int x) { return __builtin_ctz(x); }
template<>
int ctz<unsigned long>(unsigned long x) { return __builtin_ctzl(x); }
template<>
int ctz<unsigned long long>(unsigned long long x) { return __builtin_ctzll(x); }
