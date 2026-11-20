#include <algorithm>
#include <cmath>
#include <limits>

namespace {

constexpr double EPS = std::numeric_limits<double>::epsilon(); ///< underflow accuracy

/*
 * Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104,
 * Appendix A, Eq.(A3)
 */
double omega(double a, double b) noexcept
{
   const double d = a - b;
   const double C = 0.5 * (a + b) - 0.25 * d * d - 0.25;
   const double c = 0.5 * (a + b - 1);

   if (C > 0) {
      const double r = std::sqrt(C);
      return r * std::atan2(r, c);
   } else if (C < 0) {
      const double sC = std::sqrt(-C);
      return 0.5 * sC * std::log((c - sC) / (c + sC));
   }

   return 0;
}

// H(z) = z log(z) + (1-z) log(1-z).
double H(double z) noexcept
{
   if (z == 0 || z == 1) {
      return 0;
   }
   if (z > 0.5) {
      z = 1 - z;
   }
   return z * std::log(z) + (1 - z) * std::log1p(-z);
}

double log_ratio(double a, double b) noexcept
{
   const double r = a / b;

   if (r >= std::numeric_limits<double>::min() && std::isfinite(r)) {
      return std::log(r);
   }

   return std::log(a) - std::log(b);
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
double reb0xx(double s, double x, double q) noexcept
{
   const double eps_q = EPS * q;

   if (s < eps_q && x < eps_q) {
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
         + d*1./360360))))));
   } else if (s <= 4 * x) {
      const double xs = 4 * x / s;
      return 2 - std::log(x / q) -
             2 * std::sqrt(xs - 1) * std::asin(std::sqrt(1 / xs));
   } else if (s < 1e2 * x) {
      const double xs = 4 * x / s;
      const double sq = std::sqrt(1 - xs);
      return 2 - std::log(x / q) - sq * std::log((1 + sq) / (1 - sq));
   } else if (s*EPS < x) {
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

/// returns x * log(x)
double xlogx(double x) noexcept
{
   return x == 0 ? 0 : x * std::log(x);
}

} // anonymous namespace

namespace pv {

/*
 * Real part of B0 1-loop function.
 *
 * @param s squared momentum
 * @param x squared mass
 * @param y squared mass
 * @param q squared renormalization scale, q > 0
 * @return Re[B0(s, x, y, q)]
 */
double reb0(double s, double x, double y, double q) noexcept
{
   const double eps_q = EPS * q;

   // protect against infrared divergence
   if (s < eps_q && x < eps_q && y < eps_q) {
      return 0;
   }

   if (x > y) {
      std::swap(x, y);
   }

   // x == y
   if (x*(1 + EPS) > y) {
      return reb0xx(s, x, q);
   }

   if (s >= 0 && y > 0) {
      const double r = x / y;
      const double t = s / y;

      // In diesem Bereich ist A*A - 4*t > 0 und beta < etwa 0.52.
      if (r <= 1e-2 && t <= 0.5) {
         const double A = 1 - r + t;
         const double beta = 2*t / (A + std::sqrt(A*A - 4*t));
         const double delta = r / (1 - beta);

         double fbeta = 0;
         double power = beta;
         for (int n = 1; n < 128; ++n) {
            const double next = fbeta + power / (n * (n + 1.0));
            if (next == fbeta) break;
            fbeta = next;
            power *= beta;
         }

         const double mass_correction =
            delta == 0 ? 0 : delta * std::log(delta) / (1 - delta);

         return 1 - std::log(y / q) + mass_correction + fbeta;
      }
   }

   // s << y
   if (s >= 0 && y > 0 && s <= 1e-5*y) {
      const double t = s / y;
      const double r = x / y;
      const double h = 1 - r;
      double f0, f1, f2;

      if (h < 0.1) {
         f0 = f1 = f2 = 0;
         double p = 1;
         for (int k = 0; k < 16; ++k) {
            f1 += p / ((k + 2.0) * (k + 3.0));
            f2 += p * (k + 1.0)
               / ((k + 3.0) * (k + 4.0) * (k + 5.0));
            p *= h;
         }
         p = h;
         for (int n = 1; n <= 16; ++n) {
            f0 += p / (n * (n + 1.0));
            p *= h;
         }
      } else {
         const double u = xlogx(r);
         const double h2 = h * h;
         f0 = 1 + u / h;
         f1 = (1 + r) / (2 * h2) + u / (h2 * h);
         f2 = (1 + 9*r - 9*r*r - r*r*r + 6*(1+r)*u)
            / (6 * h2 * h2 * h);
      }

      return -std::log(y / q) + f0 + t * (f1 + t * f2);
   }

   if (y < 1e-4 * s) { // x <= y << s
      const double a = x / s;
      const double b = y / s;
      const double k = 1 - a - b;
      const double d = std::sqrt(k*k - 4*a*b);
      const double A = 1 + b - a;
      const double z_small = 2*b / (A + d);
      const double one_minus_z_large = 2*a / (1 - b + a + d);
      // H(z_large) == H(1 - z_large).
      return 2 - log_ratio(s, q) - H(z_small) - H(one_minus_z_large);
   } else if (x < EPS * y) { // x << y
      // s == y
      if (std::abs(y - s) < EPS * std::max(s, y)) {
         const double pi = 3.1415926535897932;
         const double a = x / s;
         return 2 - std::log(y / q) - pi * std::sqrt(a) + a - 0.5 * xlogx(a);
      } else if (s < y) {
         const double b = y / s;
         return 2 - std::log(y / q) + (b - 1) * std::log1p(-1 / b);
      } else { // s > y
         const double b = y / s;
         return 2 - std::log(y / q) + (b - 1) * std::log(1 / b - 1);
      }
   }

   const double a = x / s;
   const double b = y / s;
   const double delta = a - b;
   const double c = 0.5 * (a + b - 1);
   const double C = 0.5 * (a + b) - 0.25 * delta * delta - 0.25;

   if (C < 0 && c < 0) {
      const double w = std::sqrt(-C);
      return 2 - std::log(s / q)
         + (w - 0.5 * (1 + delta)) * std::log(a)
         + (w - 0.5 * (1 - delta)) * std::log(b)
         - 2 * w * std::log(w - c);
   }

   return 2 - std::log(s / q)
          - 0.5 * (1 + delta) * std::log(a)
          - 0.5 * (1 - delta) * std::log(b)
          - 2 * omega(a, b);
}

} // namespace pv
