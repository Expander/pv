#include "pv.hpp"
#include <cstdio>
#include <string>

int main(int argc, char* argv[])
{
   if (argc == 2) {
      const double x = std::stod(argv[1]);
      const double q = 1;
      std::printf("Re[a0(x = %.17g, q = %.17g)] = %.17g\n", x, q, pv::a0(x, q));
   } else if (argc == 4) {
      const double s = std::stod(argv[1]);
      const double x = std::stod(argv[2]);
      const double y = std::stod(argv[3]);
      const double q = 1;
      std::printf("Re[b0(s = %.17g, x = %.17g, y = %.17g, q = %.17g)] = %.17g\n", s, x, y, q, pv::b0(s, x, y, q));
   } else {
      std::printf(
         "Usage: (q = 1\n"
         "  <program-name> x\n"
         "  <program-name> s x y\n"
         );
   }

   return 0;
}
