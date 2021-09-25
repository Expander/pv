#include <cmath>

namespace {

/*
 * Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104,
 * Appendix A, Eq.(A3)
 */
double omega(double a, double b) noexcept
{
   const double C = 0.5 * (a + b) - 0.25 * (a - b)*(a - b) - 0.25;
   double res = 0;

   if (C > 0) {
      const double sC = std::sqrt(C);
      res = sC * (std::atan((1 + a - b) / (2 * sC)) +
                  std::atan((1 - a + b) / (2 * sC)));
   } else if (C < 0) {
      const double sC = std::sqrt(-C);
      res = 0.5 * sC *
            std::log((0.5 * (a + b - 1) - sC) /
                     (0.5 * (a + b - 1) + sC));
   } else {
      res = 0;
   }

   return res;
}

} // anonymous namespace

namespace pv {

/*
 * Real part of B0 1-loop function.
 *
 * @param s squared momentum
 * @param x squared mass
 * @param y squared mass
 * @param q squared renormalization scale
 * @return real part of B0
 */
double b0(double s, double x, double y, double q) noexcept
{
   double res = 0;

   if (s == 0) {
      if (x == 0 && y != 0) {
         res = 1 - std::log(y / q);
      } else if (x != 0 && y == 0) {
         res = 1 - std::log(x / q);
      } else if (x == y) {
         res = -std::log(x / q);
      } else {
         res = 1 - std::log(y / q) + x / (x - y) * std::log(y / x);
      }
   } else {
      if (x == 0 && y != 0) {
         if (y != s) {
            res = -std::log(y / q) + 2
                  + (y / s - 1) * std::log(fabs(1 - s / y));
         } else {
            res = -std::log(y / q) + 2;
         }
      } else if (y == 0 && x != 0) {
         if (x != s) {
            res = -std::log(x / q) + 2
                  + (x / s - 1) * std::log(fabs(1 - s / x));
         } else {
            res = -std::log(x / q) + 2;
         }
      } else if (y == 0 && x == 0) {
         res = -std::log(s / q) + 2;
      } else if (x == y) {
         if (s <= 4 * y) {
            res = -std::log(y / q) + 2
                  - 2 * std::sqrt(4 * y / s - 1)
                      * std::asin(std::sqrt(s / (4 * y)));
         } else {
            const double sC = std::sqrt(1 - 4 * y / s);
            res = -std::log(y / q) + 2
                  + sC * std::log(s * (1 - sC) / (2 * y) - 1);
         }
      } else {
         const double a = x / s;
         const double b = y / s;
         const double delta = a - b;
         res = -std::log(s / q) + 2
               - 0.5 * (1 + delta) * std::log(a)
               - 0.5 * (1 - delta) * std::log(b)
               - 2 * omega(a, b);
      }
   }

   return res;
}

} // namespace pv
