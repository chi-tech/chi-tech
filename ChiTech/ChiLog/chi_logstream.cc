#include "chi_logstream.h"

LogStream::~LogStream()
{
  if (dummy) return;
//  std::cout << "Sa" << std::flush;
  std::string line, oline;
  while (std::getline(*this, line))
    oline += log_header + line + '\n';

  if (!oline.empty())
    *log_stream << oline << std::flush;

//  std::cout << "Ba" << std::flush;
}