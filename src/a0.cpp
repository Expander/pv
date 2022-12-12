#include <cmath>

namespace pv {

/*
 * Real part of A0 1-loop function.
 *
 * @param x squared mass
 * @param q squared renormalization scale
 * @return real part of A0
 */
double a0(double x, double q) noexcept
{
   if (q == 0) {
      return 0;
   }

   if (x == 0) {
      return 0;
   }
   return x * (1 - std::log(x/q));
}

} // namespace pv
