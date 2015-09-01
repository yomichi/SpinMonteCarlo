#ifndef DETAIL_LOGGING_MACRO_HPP
#define DETAIL_LOGGING_MACRO_HPP

#ifdef ENABLE_LOGGING

#define NV(name) (#name) << " = " << (name) << " "

#define LOGGING(logger, args) do{util::WriterPtr(new util::Writer(logger)) << args; }while(false)

#define SCOPED_LOGGING(logger, s) std::stringstream ss_##__LINE__; ss_##__LINE__ << s; util::ScopedLogger scopedlogger_##__LINE__(logger, ss_##__LINE__.str())

#else

#define LOGGING(logger, args) (0)
#define SCOPED_LOGGING(logger, str) (0)
#endif

#ifdef ENABLE_LOGGING_ASSERTION
#ifdef ENABLE_LOGGING
#define logging_assert(logger, x) \
  do{\
    if(!(x)){\
      LOGGING(logger, \
          "assertion failed at " << __FILE__ << " (" << __LINE__ << "): " << #x \
          );\
      std::cerr << "assertion failed at " << __FILE__ << " (" << __LINE__ << "): " << #x << std::endl; \
      std::exit(127);\
    }\
  }while(false)
#else // ENABLE_LOGGING
#define logging_assert(logger, x) \
  do{\
    if(!(x)){\
      std::cerr << "assertion failed at " << __FILE__ << " (" << __LINE__ << "): " << #x << std::endl; \
      std::exit(127);\
    }\
  }while(false)
#endif // ENABLE_LOGGING
#else // ENABLE_LOGGING_ASSERTION
#define logging_assert(logger, x) (0)
#endif // ENABLE_LOGGING_ASSERTION

#endif // DETAIL_LOGGING_MACRO_HPP
