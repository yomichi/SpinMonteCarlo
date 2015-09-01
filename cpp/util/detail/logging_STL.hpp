#ifndef DETAIL_LOGGING_STL_HPP
#define DETAIL_LOGGING_STL_HPP

#ifdef ENABLE_LOGGING

#include <vector>
#include <list>
#include <set>

namespace util{


#define STL_CONTAINER_ONE(name) \
template <class T, class Allocator> \
inline WriterPtr operator<<(WriterPtr writer, name<T, Allocator> const& xs) \
{ \
  if(xs.empty()){\
    writer << "[]";\
    return writer; \
  } \
  writer->print("["); \
  typename name<T, Allocator>::const_iterator iter = xs.begin(); \
  typename name<T, Allocator>::const_iterator end = xs.end(); \
  writer << (*iter); \
  ++iter; \
  for(;iter != end; ++iter){ \
    writer << " " << (*iter); \
  }\
  writer << "]"; \
  return writer;\
}

STL_CONTAINER_ONE(std::vector);
STL_CONTAINER_ONE(std::list);
STL_CONTAINER_ONE(std::set);

#undef STL_CONTAINER_ONE

} // end of namespace util

#endif // ENABLE_LOGGING
#endif // DETAIL_LOGGING_STL_HPP
