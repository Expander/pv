#pragma once
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

const char PATH_SEPARATOR =
#ifdef _WIN32
   '\\';
#else
   '/';
#endif

/**
 * Reads real numbers from a file line by line into a vector of vectors of T.
 *
 * @param filename file name
 * @tparam T data type
 *
 * @return vector of vectors of doubles.
 */
std::vector<std::vector<double>>
read_from_file(const std::string& filename)
{
   std::vector<std::vector<double>> data;
   std::string line;
   std::ifstream fstr(filename);

   if (!fstr.is_open()) {
      throw std::runtime_error("Cannot open file: " + filename);
   }

   while (std::getline(fstr, line)) {
      std::istringstream iss(line);
      std::vector<double> vec{std::istream_iterator<double>(iss),
                              std::istream_iterator<double>()};
      data.emplace_back(vec);
   }

   return data;
}
