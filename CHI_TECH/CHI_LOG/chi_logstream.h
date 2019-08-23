#ifndef _chi_logstream_h
#define _chi_logstream_h

#include <iostream>
#include <sstream>

//###################################################################
/** Log stream for adding header information to a string stream.*/
class CHI_LOG_STREAM : public std::stringstream
{
private:
  std::ostream* log_stream;
  std::string  log_header;

public:
  /** Creates a string stream.*/
  CHI_LOG_STREAM(std::ostream* output_stream, std::string header)
  : log_stream(output_stream),
    log_header(std::move(header))
  { }

  /** Flushes the broken-up/headered stream to the output.*/
  virtual ~CHI_LOG_STREAM()
  {
    std::string line, oline;
    while (std::getline(*this, line))
    {
      oline += log_header + line + '\n';
    }

    if (!oline.empty())
      *log_stream << oline << std::flush;
  }

  CHI_LOG_STREAM(const CHI_LOG_STREAM& other)
  {
    log_stream = other.log_stream;
    log_header = other.log_header;
  }
};

struct DummyStream: public std::ostream
{
  struct DummyStreamBuffer : std::streambuf
  {
    virtual int overflow(int c) { return c; };
  } buffer;

  DummyStream(): std::ostream(&buffer) {}
};

#endif