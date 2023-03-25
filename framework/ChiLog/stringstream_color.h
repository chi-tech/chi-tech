#ifndef CHITECH_STRINGSTREAM_COLOR_H
#define CHITECH_STRINGSTREAM_COLOR_H

#include <string>

namespace chi_objects
{
enum StringSteamColorCode
{
  RESET = 0,
  FG_BOLD = 1,
  FG_UNDERLINE = 4,
  FG_BOLD_OFF = 21,
  FG_UNDERLINE_OFF = 24,
  FG_RED = 31,
  FG_GREEN = 32,
  FG_YELLOW = 33,
  FG_BLUE = 34,
  FG_MAGENTA = 35,
  FG_CYAN = 36,
  FG_WHITE = 37,
  FG_DEFAULT = 39,
};

inline std::string StringStreamColor(StringSteamColorCode code)
{
  return std::string("\033[") + std::to_string(code) + "m";
}

}//chi_objects

#endif //CHITECH_STRINGSTREAM_COLOR_H
