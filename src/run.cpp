#include "pv.hpp"
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <sstream>
#include <string>


double string_to_double(char* str)
{
   std::istringstream iss(str);
   double x;
   iss >> x;
   if (iss.fail()) {
      throw std::range_error(std::string("input not convertible to double: ") + str);
   }
   return x;
}


int main(int argc, char* argv[])
{
   if (argc == 2) {
      // A0 function
      const double x = string_to_double(argv[1]);
      const double q = 1;
      std::printf("Re[a0(x = %.17g, q = %.17g)] = %.17g\n",
                  x, q, pv::rea0(x, q));
   } else if (argc == 4) {
      // B0 function
      const double s = string_to_double(argv[1]);
      const double x = string_to_double(argv[2]);
      const double y = string_to_double(argv[3]);
      const double q = 1;
      std::printf("Re[b0(s = %.17g, x = %.17g, y = %.17g, q = %.17g)] = %.17g\n",
                  s, x, y, q, pv::reb0(s, x, y, q));
   } else {
      std::printf(
         "Usage: <program-name> [parameters]\n"
         "\n"
         "This program calculates n-point 1-loop integrals.\n"
         "The value of n depends on the number of given parameters:\n"
         "\n"
         "  n = 1: <program-name> m^2\n"
         "  n = 2: <program-name> p^2 m1^2 m2^2\n"
         "\n"
         "The parameters are always interpreted as squared masses or squared energies.\n"
         "The renormalization scale Q is always set to Q = 1.\n"
         );
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
