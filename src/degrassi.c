#include <math.h>


/*
 * Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104,
 * Appendix A, Eq.(A3)
 */
static double omega(double a, double b)
{
   const double C = 0.5 * (a + b) - 0.25 * (a - b)*(a - b) - 0.25;
   double res = 0;

   if (C > 0) {
      const double sC = sqrt(C);
      res = sC * (atan((1 + a - b) / (2 * sC)) +
                  atan((1 - a + b) / (2 * sC)));
   } else if (C < 0) {
      const double sC = sqrt(-C);
      res = 0.5 * sC *
            log((0.5 * (a + b - 1) - sC) /
                (0.5 * (a + b - 1) + sC));
   } else {
      res = 0;
   }

   return res;
}


/*
 * Real part of B0 1-loop function.
 *
 * From Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104,
 * Appendix A, Eq.(A2)
 *
 * @param p2 squared momentum
 * @param m1 squared mass
 * @param m2 squared mass
 * @param q2 squared renormalization scale
 * @return real part of B0
 */
double b0_degrassi(double p2, double m1, double m2, double q2)
{
   double res = 0;

   if (p2 == 0) {
      if (m1 == 0 && m2 != 0) {
         res = 1 - log(m2 / q2);
      } else if (m1 != 0 && m2 == 0) {
         res = 1 - log(m1 / q2);
      } else if (m1 == m2) {
         res = -log(m1 / q2);
      } else {
         res = 1 - log(m2 / q2) + m1 / (m1 - m2) * log(m2 / m1);
      }
   } else {
      if (m1 == 0 && m2 != 0) {
         if (m2 != p2) {
            res = -log(m2 / q2) + 2 + (m2 / p2 - 1) * log(fabs(1 - p2 / m2));
         } else {
            res = -log(m2 / q2) + 2;
         }
      } else if (m2 == 0 && m1 != 0) {
         if (m1 != p2) {
            res = -log(m1 / q2) + 2 + (m1 / p2 - 1) * log(fabs(1 - p2 / m1));
         } else {
            res = -log(m1 / q2) + 2;
         }
      } else if (m2 == 0 && m1 == 0) {
         res = -log(p2 / q2) + 2;
      } else {
         const double delta = (m1 - m2)/p2;
         const double a = m1 / p2;
         const double b = m2 / p2;
         res = -log(p2 / q2) + 2
               - 0.5 * (1 + delta) * log(a)
               - 0.5 * (1 - delta) * log(b)
               - 2 * omega(a, b);
      }
   }

   return res;
}
