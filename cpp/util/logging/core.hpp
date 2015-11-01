#ifndef LOGGING_CORE_HPP
#define LOGGING_CORE_HPP

#include <string>
#include <sstream>

#ifdef ENABLE_LOGGING

#include <iostream>
#include <fstream>
#include <boost/shared_ptr.hpp>

namespace util{
class Writer;
class Logger{
private:
  std::string filename;
  boost::shared_ptr<std::ofstream> ofs;
public:
  explicit Logger(const char* logger_id)
  {
    init(std::string(logger_id));
  }
  explicit Logger(std::string const& logger_id)
  {
    init(logger_id);
  }
  void init(std::string const& logger_id)
  {
    filename = logger_id + ".log";
    ofs = boost::shared_ptr<std::ofstream>(new std::ofstream(filename.c_str()));
  }
  void reset(){ ofs->close(); ofs->open(filename.c_str());}
  friend class Writer;
};

class Writer{
private:
  boost::shared_ptr<std::ofstream> ofs;
public:
  Writer(Logger &logger):ofs(logger.ofs){}
  ~Writer(){ (*ofs) << std::endl;}
  template <class T>
  void print(T const& x){ (*ofs) << x;}
  void manipulate(std::ostream& (*manip)(std::ostream&)){ (*ofs) << manip;}
};

typedef boost::shared_ptr<Writer> WriterPtr;

template <class T>
inline WriterPtr operator<<(WriterPtr writer, T const& x)
{
  writer->print(x);
  return writer;
}
inline WriterPtr operator<<(WriterPtr writer, std::ostream& (*manip)(std::ostream&))
{
  writer->manipulate(manip);
  return writer;
}

template <class T>
inline WriterPtr operator<<(Logger &logger, T const& x)
{
  WriterPtr writer(new Writer(logger));
  writer << x;
  return writer;
}


class ScopedLogger {
private:
  Logger logger_;
  std::string str_;
  ScopedLogger(ScopedLogger const&);
  ScopedLogger& operator=(ScopedLogger const&);
public:
  ScopedLogger(Logger const &logger, std::string const& str):logger_(logger), str_(str)
  {
    LOGGING(logger_, str_ << " in  ");
  }
  ~ScopedLogger()
  {
    LOGGING(logger_, str_ << " out ");
  }
};
} // end of namespace util

#else // ENABLE_LOGGING
namespace util{
class Logger{
public:
  explicit Logger(const char*){}
  explicit Logger(std::string const&){}
  void reset(){}
};
} // end of namespace util
#endif // ENABLE_LOGGING

#endif // LOGGING_CORE_HPP
