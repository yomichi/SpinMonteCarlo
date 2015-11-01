#ifndef LOGGING_MACRO_HPP
#define LOGGING_MACRO_HPP

#ifdef ENABLE_LOGGING

#define NV(name) (#name) << " = " << (name) << " "

#define LOGGING(logger, args) do{util::WriterPtr(new util::Writer(logger)) << args; }while(false)

#define SCOPED_LOGGING__(logger, s,l) std::stringstream ss_##l; ss_##l << s; util::ScopedLogger scopedlogger_##l(logger, ss_##l.str())
#define SCOPED_LOGGING_(logger, s, l) SCOPED_LOGGING__(logger, s, l)
#define SCOPED_LOGGING(logger, s) SCOPED_LOGGING_(logger, s, __LINE__)

#ifdef ENABLE_LOGGING_VERBOSE
#define VERBOSE(logger, args) do{util::WriterPtr(new util::Writer(logger)) << args; }while(false)
#define SCOPED_VERBOSE(logger, s) SCOPED_LOGGING_(logger, s, __LINE__)
#else 
#define VERBOSE(logger, args) (0)
#define SCOPED_VERBOSE(logger, str) (0)
#endif // ENABLE_LOGGING_VERBOSE

#else // ENABLE_LOGGING

#define LOGGING(logger, args) (0)
#define VERBOSE(logger, args) (0)
#define SCOPED_LOGGING(logger, str) (0)
#define SCOPED_VERBOSE(logger, str) (0)
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

#endif // LOGGING_MACRO_HPP
