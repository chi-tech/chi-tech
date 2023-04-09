#include "stringstream_color.h"

#include "chi_runtime.h"

std::string chi_objects::StringStreamColor(StringSteamColorCode code)
{
  if (chi::run_time::suppress_color_) return {};
  return std::string("\033[") + std::to_string(code) + "m";
}