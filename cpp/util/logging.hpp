#ifndef LOGGING_HPP
#define LOGGING_HPP

/*
 * usage
 *

 To switch on logger, define ENABLE_LOGGING
 
 // initialize logger
 // log messages will be written in './<logger_id>.log'.
 util::Logger logger("logger_id");

 LOGGING(logger, "log"); // "log"

 int x = 42;
 double y = 3.14;
 LOGGING(logger, NV(x) << NV(y)); // "x = 42 y = 3.14"

 // clear log file
 logger.reset();

 // logger can be copied
 util::Logger logger2 = logger;

 {
   SCOPED_LOGGING(logger, "scope"); // "scope in"
 } // "scope out"

 // When ENABLE_LOGGING_ASSERTION is defined,
 // the following assertion message
 // "assertion failed at <filename> (<line number>): 1 == 0"
 // and the program will stop soon.
 logging_assert(logger, 1 == 0);


 // VERBOSE is the almost same as LOGGING,
 // but this prints only when ENABLE_LOGGING_VERBOSE is defined
 VERBOSE(logger, "verbose output");
 { SCOPED_VERBOSE(logger, "verbose output"); }

 * end of usage
 */

#include "logging/macro.hpp"
#include "logging/core.hpp"
#include "logging/STL.hpp"

#endif // LOGGING_HPP
