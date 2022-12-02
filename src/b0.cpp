#include <cmath>
#include <algorithm>

namespace {

constexpr double EPSTOL = 1.0e-15; ///< underflow accuracy

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

/**
 * Re(B0(s,x,x,q)), Eq.(2.4) from [hep-ph/0701051]
 *
 * @param s squared momentum
 * @param x squared mass
 * @param q squared renormalization scale
 *
 * @return Re(B0(s,x,x,q))
 */
double b0xx(double s, double x, double q) noexcept
{
   if (s < EPSTOL * q && x < EPSTOL * q) {
      return 0; // IR divergence
   } else if (s < 1e-3*x) {
      const double d = s / x;
      return -std::log(x / q)
         + d*(1./6
         + d*(1./60
         + d*(1./420
         + d*(1./2520
         + d*(1./13860
         + d*(1./72072
         + d/360360))))));
   } else if (s <= 4 * x) {
      return 2 - std::log(x / q) -
             2 * std::sqrt(4 * x / s - 1) *
                std::asin(std::sqrt(s / (4 * x)));
   } else if (s < 1e2 * x) {
      const double sq = std::sqrt(1 - 4 * x / s);
      return 2 - std::log(x / q) +
         sq * std::log(s * (1 - sq) / (2 * x) - 1);
   } else if (s*EPSTOL < x) {
      const double d = x / s;
      const double logd = std::log(d);
      return 2 - std::log(s / q)
         + d * (2 * (1 - logd)
         + d * (-1 - 2 * logd
         + d * (-10./3 - 4 * logd
         + d * (-59./6 - 10 * logd
         + d * (-449./15 - 28 * logd
         + d * (-1417./15 - 84 * logd
         + d * (-32254./105 - 264 * logd
         + d * (-429697./420 - 858 * logd))))))));
   } else {
      return 2 - std::log(s / q);
   }
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
   const double eps_q = EPSTOL * q;

   // protect against infrared divergence
   if (s < eps_q && x < eps_q && y < eps_q) {
      return 0;
   }

   if (x > y) {
      std::swap(x, y);
   }

   // x == y
   if (x*(1 + EPSTOL) > y) {
      return b0xx(s, x, q);
   }

   // s << x
   if (s <= 1e-11 * x) {
      if (x < EPSTOL * y) {
         return 1 - std::log(y / q);
      }

      return 1 - std::log(y / q) + x * std::log(y / x) / (x - y);
   }

   double res = 0;

   if (y < 1e-4 * s) { // x < y << s
      res = -std::log(s / q) + 2;
   } else if (x < 1e-3 * y) { // x << y
      // s == y
      if (std::abs(y - s) < EPSTOL * std::max(s, y)) {
         res = -std::log(y / q) + 2;
      } else {
         res = -std::log(y / q) + 2
               + (y / s - 1) * std::log(std::abs(1 - s / y));
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

   return res;
}

} // namespace pv
