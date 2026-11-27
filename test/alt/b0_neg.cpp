#include <algorithm>
#include <cmath>

namespace {

/// returns x*log(|x|)
double xlogabsx(double x) noexcept
{
   return x == 0 ? 0 : x*std::log(std::abs(x));
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
   return x*std::log(std::abs(x)) + (1 - x)*std::log(std::abs(1 - x));
}

double lambda(double s, double x, double y) noexcept
{
   const double l = s - x - y;
   return l*l - 4*x*y;
}

double reb0_x0(double s, double y, double q) noexcept
{
   if (y == 0)
      return s == 0 ? 0 : 2 - std::log(std::abs(s/q));

   if (s == 0)
      return 1 - std::log(std::abs(y/q));

   if (s == y)
      return 2 - std::log(std::abs(y/q));

   // |s/y| klein: konvergente Reihe, auch für negatives s/y.
   if (std::abs(s) <= 0.5*std::abs(y)) {
      const double t = s/y;
      double sum = 0;
      double power = t;
      for (int n = 1; n < 128; ++n) {
         const double next = sum + power/(n*(n + 1.0));
         if (next == sum) break;
         sum = next;
         power *= t;
      }
      return 1 - std::log(std::abs(y/q)) + sum;
   }

   // |y/s| klein: Darstellung über H_R(y/s).
   if (std::abs(y) <= 0.5*std::abs(s)) {
      const double b = y/s;
      const double h = xlogabsx(b) + (1 - b)*std::log1p(-b);
      return 2 - std::log(std::abs(s/q)) - h;
   }

   // Vergleichbare Beträge: nahe s == y die Differenz direkt bilden.
   if (std::signbit(s) == std::signbit(y)) {
      const double d = y - s;
      return 2 - std::log(std::abs(y/q))
             + (d/s)*std::log(std::abs(d/y));
   }

   // Entgegengesetzte Vorzeichen: hier gibt es keine Auslöschung in 1-s/y.
   return 2 - std::log(std::abs(y/q))
          + (y/s - 1)*std::log1p(-s/y);
}

} // anonymous namespace

namespace pv {

double b0_neg(double s, double x, double y, double q) noexcept
{
   const double lam = lambda(s, x, y);

   if (std::abs(x) > std::abs(y)) {
      std::swap(x, y);
   }

   if (x == 0) {
      return reb0_x0(s, y, q);
   }

   if (s != 0) {
      const double u = (s + y - x)/(2*s);
      const double v = std::sqrt(std::abs(lam))/(2*std::abs(s));

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
