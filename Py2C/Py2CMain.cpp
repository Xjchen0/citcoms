#include "Py2C.hpp"
#include <filesystem>
#include <string_view>

namespace fs = std::filesystem;

constexpr std::string_view usage = 
  "Usage: Py2C infile outfile [save_all=False]\n\n"
  "Where:\n"
  "  infile        Name of Python/INI configuration file (.cfg)\n"
  "  outfile       Name of C/C++ version of the configuration file\n"
  "  save_all      Optional. 'True' writes all known CitcomS parameters grouped by section.\n"
  "                'False' (default) writes only parameters defined in infile.\n";

int main(int argc, char* argv[])
{
  if (argc < 3 || argc > 4)
  {
    std::cerr << usage << std::endl;
    return 1;
  }

  const fs::path inpath(argv[1]);
  if (!fs::exists(inpath))
  {
    std::cerr << "Error: Input file does not exist: " << inpath << std::endl;
    return 1;
  }

  Py2CConverter py2c(argv[1], argv[2]);

  bool saveall = false;
  if (argc == 4)
  {
    std::string saveallstr(argv[3]);
    std::transform(saveallstr.begin(), saveallstr.end(), saveallstr.begin(), [](unsigned char c) {
      return std::tolower(c);
    });
    saveall = (saveallstr == "true" || saveallstr == "1" || saveallstr == "yes");
  }

  py2c.convert(saveall);
  return 0;
}
