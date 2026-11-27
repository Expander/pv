#include <cmath>

namespace {

/// returns x*log(x)
double xlogx(double x) noexcept
{
   return x == 0 ? 0 : x*std::log(x);
}

/*
 * H(x) = x*log(x) + (1-x)*log(1-x).
 * H(z) == H(1 - z).
 */
double H(double x) noexcept
{
   if (x == 0 || x == 1) {
      return 0;
   }
   if (x > 0.5) {
      x = 1 - x;
   }
   return x*std::log(x) + (1 - x)*std::log1p(-x);
}

double lambda(double s, double x, double y) noexcept
{
   const double l = s - x - y;
   return l*l - 4*x*y;
}

} // anonymous namespace

namespace pv {

double b0_neg(double s, double x, double y, double q) noexcept
{
   const double lam = lambda(s, x, y);
   const double u = (s + y - x)/(2*s);
   const double v = std::sqrt(std::abs(lam))/(2*s);

   if (s != 0) {
      if (lam >= 0) {
         return 2 - std::log(std::abs(s/q)) - H(u + v) - H(u - v);
      } else { // lam < 0
         return 2 - std::log(std::abs(s/q))
            - 2*(u*std::log(std::hypot(u,v))
                 + (1 - u)*std::log(std::hypot(1 - u,v))
                 + v*(std::atan2(u,v) + std::atan2(1 - u,v)));
      }
   } else { // s == 0
      if (x != y) {
         return 1 - (x*std::log(std::abs(x/q)) - y*std::log(std::abs(y/q)))/(x - y);
      } else if (y != 0) { // x == y and y != 0
         return -std::log(std::abs(x/q));
      } else { // x == y == 0
         return 0;
      }
   }
}

} // namespace pv
