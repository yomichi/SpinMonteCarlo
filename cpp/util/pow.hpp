#ifndef UTIL_POW_HPP
#define UTIL_POW_HPP

namespace util{

template <class T>
inline T pow2(const T &x) { return x * x;}

template <class T>
inline T pow3(const T &x) { return x * x * x;}

template <class T>
inline T pow4(const T &x) { return pow2(pow2(x));}

template <class T>
inline T pow5(const T &x) { return pow4(x) * x;}

template <class T>
inline T pow6(const T &x) { return pow2(pow3(x));}

} // end of namespace util

#endif
