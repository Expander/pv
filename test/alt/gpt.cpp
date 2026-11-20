#include <algorithm>
#include <cmath>
#include <limits>

namespace {

constexpr double EPS = std::numeric_limits<double>::epsilon();

double log_ratio(double a, double b) noexcept
{
   // a > 0, b > 0. Vermeidet Über- oder Unterlauf bei a/b.
   const double r = a / b;
   if (r >= std::numeric_limits<double>::min() && std::isfinite(r))
      return std::log(r);

   return std::log(a) - std::log(b);
}

double xlogx(double x) noexcept
{
   return x == 0 ? 0 : x * std::log(x);
}

/*
 * G(v) = sum_{n>=1} v^n/[n(n+1)]
 *      = 1 + (1-v)/v * log(1-v),  0 <= v <= 1.
 *
 * distance_to_one wird gegebenenfalls separat und genauer berechnet,
 * statt 1-v durch Auslöschung zu bilden.
 */
double G(double v, double distance_to_one) noexcept
{
   if (v == 0)
      return 0;

   if (v <= 0.5) {
      double sum = 0;
      double power = v;

      for (int n = 1; n < 256; ++n) {
         const double next = sum + power / (n * (n + 1.0));
         if (next == sum)
            break;

         sum = next;
         power *= v;
      }
      return sum;
   }

   if (distance_to_one == 0)
      return 1;

   return 1 + distance_to_one / v * std::log(distance_to_one);
}

/*
 * H(z) = z log(z) + (1-z) log(1-z).
 * Durch log1p auch für kleines z stabil.
 */
double H(double z) noexcept
{
   if (z == 0 || z == 1)
      return 0;

   if (z > 0.5)
      z = 1 - z;

   return xlogx(z) + (1 - z) * std::log1p(-z);
}

/*
 * B0 für x = 0 und y > 0; s > 0.
 */
double reb0_massless_x(double s, double y, double q) noexcept
{
   if (s <= y) {
      const double t = s / y;
      const double distance_to_one = (y - s) / y;
      return 1 - log_ratio(y, q) + G(t, distance_to_one);
   }

   const double b = y / s;

   // H(b) = H(1-b). Bei b nahe 1 die kleine Differenz
   // aus s-y statt aus 1-b berechnen.
   const double small_argument =
      b > 0.5 ? (s - y) / s : b;

   return 2 - log_ratio(s, q) - H(small_argument);
}

/*
 * Direkte Integration des glatten Feynman-Parameter-Integranden
 * mit einer 16-Punkt-Gauß-Legendre-Regel.
 *
 * Aufruf nur für x/y >= 0.5 und s/y <= 0.1.
 */
double reb0_small_s_near_equal(double s, double x,
                              double y, double q) noexcept
{
   static const double nodes[8] = {
      0.989400934991649933,
      0.944575023073232576,
      0.865631202387831744,
      0.755404408355003034,
      0.617876244402643748,
      0.458016777657227386,
      0.281603550779258913,
      0.0950125098376374402
   };

   static const double weights[8] = {
      0.0271524594117540949,
      0.0622535239386478929,
      0.0951585116824927848,
      0.124628971255533872,
      0.149595988816576732,
      0.169156519395002538,
      0.182603415044923589,
      0.189450610455068496
   };

   const double h = (y - x) / y;
   const double t = s / y;

   double integral = 0;

   for (int i = 0; i < 8; ++i) {
      const double z1 = (1 - nodes[i]) / 2;
      const double z2 = (1 + nodes[i]) / 2;

      const double f1 =
         std::log1p(-h*z1 - t*z1*(1 - z1));
      const double f2 =
         std::log1p(-h*z2 - t*z2*(1 - z2));

      integral += weights[i] * (f1 + f2);
   }

   return -log_ratio(y, q) - integral / 2;
}

} // anonymous namespace

namespace pv {

/*
 * Realteil von B0(s,x,y,q).
 *
 * s: Impulsquadrat, s >= 0
 * x,y: Massenquadrate, x,y >= 0
 * q: Quadrat der Renormierungsskala, q > 0
 *
 * Alle Eingaben müssen endlich sein.
 */
double b0_gpt_astra(double s, double x, double y, double q) noexcept
{
   if (!std::isfinite(s) || !std::isfinite(x) ||
       !std::isfinite(y) || !std::isfinite(q) ||
       s < 0 || x < 0 || y < 0 || q <= 0) {
      return std::numeric_limits<double>::quiet_NaN();
   }

   if (x > y)
      std::swap(x, y);

   // Beide Massen verschwinden.
   if (y == 0) {
      if (s == 0)
         return std::numeric_limits<double>::infinity();

      return 2 - log_ratio(s, q);
   }

   // Exakt verschwindender Impuls: keine Schwellenberechnung nötig.
   if (s == 0) {
      const double h = (y - x) / y;
      const double r = x / y;
      return -log_ratio(y, q) + G(h, r);
   }

   if (x == 0)
      return reb0_massless_x(s, y, q);

   const double r = x / y;
   const double t = s / y;

   // Ist sqrt(x/y) unterhalb der relevanten double-Auflösung,
   // ist selbst die maximale Schwellenkorrektur nicht auflösbar.
   // sqrt(x)/sqrt(y) vermeidet vorherigen Unterlauf von x/y.
   const double sqrt_r = std::sqrt(x) / std::sqrt(y);
   if (sqrt_r <= EPS / 4)
      return reb0_massless_x(s, y, q);

   // Bei nahezu gleichen Massen und kleinem Impuls würden die
   // geschlossenen Formeln Zahlen der Größe ~2 subtrahieren,
   // um einen kleinen Impulsbeitrag zu erhalten.
   if (r >= 0.5 && t <= 0.1)
      return reb0_small_s_near_equal(s, x, y, q);

   const double lower = (1 - sqrt_r) * (1 - sqrt_r);
   const double upper = (1 + sqrt_r) * (1 + sqrt_r);

   // Bereich unterhalb der Pseudoschwelle:
   // reelle Faktorisierung in (1-alpha*z)(1-beta*z).
   if (t <= lower) {
      const double D =
         std::sqrt(std::max(0.0, (upper - t) * (lower - t)));

      const double A = 1 - r + t;
      const double B = 1 + r - t;

      const double alpha = (A + D) / 2;
      const double beta  = 2*t / (A + D);

      // 1-alpha und 1-beta ohne Auslöschung berechnen.
      const double one_minus_alpha = 2*r / (B + D);
      const double one_minus_beta  = r / one_minus_alpha;

      return -log_ratio(y, q)
         + G(alpha, one_minus_alpha)
         + G(beta, one_minus_beta);
   }

   // Zwischen Pseudo- und Zweiteilchenschwelle.
   if (t < upper) {
      const double d =
         std::sqrt(std::max(0.0, (upper - t) * (t - lower)));

      /*
       * Zusammengefasste Logarithmen:
       *
       * 2 - log(y/q)
       *   + (y-s-x)/(2s) * log(x/y)
       *   - d/t * atan2(d, 1+r-t).
       *
       * (y-s)/y - r erhält insbesondere bei s ~= y und x << y
       * den kleinen Koeffizienten vor log(x/y).
       */
      const double log_coefficient =
         ((y - s) / y - r) / (2*t);

      return 2 - log_ratio(y, q)
         + log_coefficient * log_ratio(x, y)
         - (d / t) * std::atan2(d, 1 + r - t);
   }

   // Auf und oberhalb der Zweiteilchenschwelle: zwei reelle
   // Nullstellen des Feynman-Parameter-Polynoms.
   const double a = x / s;
   const double b = y / s;

   double discriminant;
   if (t < 8) {
      // Besonders nahe an der Schwelle genauer als
      // (1-a-b)^2 - 4ab.
      discriminant =
         ((t - upper) / t) * ((t - lower) / t);
   } else {
      const double k = 1 - a - b;
      discriminant = k*k - 4*a*b;
   }

   const double d = std::sqrt(std::max(0.0, discriminant));
   const double A = 1 + b - a;

   const double z_small = 2*b / (A + d);

   // Abstand der kleineren Nullstelle von 1; ohne 1-z_small.
   const double one_minus_z_small = (1 - b + a + d) / 2;

   // Produkt (1-z_small)(1-z_large) = a.
   const double one_minus_z_large = a / one_minus_z_small;

   const double h_small =
      z_small <= 0.5 ? H(z_small) : H(one_minus_z_small);

   return 2 - log_ratio(s, q)
      - h_small
      - H(one_minus_z_large);
}

} // namespace pv
