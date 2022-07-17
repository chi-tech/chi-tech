#include "chi_logstream.h"

chi_objects::LogStream::~LogStream()
{
  if (dummy) return;

  std::string line, oline;
  while (std::getline(*this, line))
    oline += log_header + line + '\n';

  if (!oline.empty())
    *log_stream << oline << std::flush;
}