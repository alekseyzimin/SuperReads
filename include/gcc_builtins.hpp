#ifndef __GCC_BUILTINS_HPP__
#define __GCC_BUILTINS_HPP__

template<typename T>
int ctz(T x);

template<>
int ctz<unsigned int>(unsigned int x);
template<>
int ctz<unsigned long>(unsigned long x);
template<>
int ctz<unsigned long long>(unsigned long long x);


#endif /* __GCC_BUILTINS_HPP__ */
