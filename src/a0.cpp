#include <cmath>
#include <limits>

namespace pv {

/*
 * Real part of A0 1-loop function.
 *
 * @param x squared mass
 * @param q squared renormalization scale, q > 0
 * @return real part of A0
 */
double rea0(double x, double q) noexcept
{
   if (q <= 0) {
      return std::numeric_limits<double>::quiet_NaN();
   }

   if (x == 0) {
      return 0;
   }

   return x * (1 - std::log(x/q));
}

} // namespace pv
