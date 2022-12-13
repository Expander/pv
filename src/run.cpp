#include "pv.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char* argv[])
{
   if (argc == 2) {
      // A0 function
      const double x = std::stod(argv[1]);
      const double q = 1;
      std::printf("Re[a0(x = %.17g, q = %.17g)] = %.17g\n", x, q, pv::rea0(x, q));
   } else if (argc == 4) {
      // B0 function
      const double s = std::stod(argv[1]);
      const double x = std::stod(argv[2]);
      const double y = std::stod(argv[3]);
      const double q = 1;
      std::printf("Re[b0(s = %.17g, x = %.17g, y = %.17g, q = %.17g)] = %.17g\n", s, x, y, q, pv::reb0(s, x, y, q));
   } else {
      std::printf(
         "Usage:\n"
         "  <program-name> m^2\n"
         "  <program-name> p^2 m1^2 m2^2\n"
         "Note: The renormalization scale Q is always set to Q = 1.\n"
         );
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
