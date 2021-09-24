#include <cmath>
#include <complex>

namespace {

constexpr double EPSTOL = 1.0e-11; ///< underflow accuracy

/// fast implementation of complex logarithm
template <class T>
std::complex<T> fast_log(const std::complex<T>& z) noexcept
{
   const T rz = std::real(z);
   const T iz = std::imag(z);

   return std::complex<T>(0.5*std::log(rz*rz + iz*iz), std::atan2(iz, rz));
}

constexpr double dabs(double a) noexcept { return a >= 0. ? a : -a; }

// can be made constexpr in C++20
double fB(const std::complex<double>& x) noexcept
{
   const double re = std::real(x);
   const double im = std::imag(x);

   if ((std::abs(re) == 0.0 || std::abs(re) == 1.0) && im == 0.0) {
      return -1.0;
   }

   return std::real(-1.0 + fast_log(1.0 - x) - x*fast_log(1.0 - 1.0/x));
}

/// fB(xp) + fB(xm)
double fB(const std::complex<double>& xp, const std::complex<double>& xm) noexcept
{
   const double rep = std::real(xp);
   const double imp = std::imag(xp);

   if ((std::abs(rep) == 0.0 || std::abs(rep) == 1.0) && imp == 0.0) {
      return -1.0 + fB(xm);
   }

   const double rem = std::real(xm);
   const double imm = std::imag(xm);

   if ((std::abs(rem) == 0.0 || std::abs(rem) == 1.0) && imm == 0.0) {
      return -1.0 + fB(xp);
   }

   return std::real(-2.0 + fast_log((1.0 - xp)*(1.0 - xm))
      - xp*fast_log(1.0 - 1.0/xp) - xm*fast_log(1.0 - 1.0/xm));
}

constexpr bool is_close(double m1, double m2, double tol) noexcept
{
   const double mmax = std::max(dabs(m1), dabs(m2));
   const double mmin = std::min(dabs(m1), dabs(m2));
   const double max_tol = tol * mmax;

   if (max_tol == 0.0 && mmax != 0.0 && tol != 0.0)
      return mmax - mmin <= tol;

   return mmax - mmin <= max_tol;
}

constexpr bool is_zero(double m, double tol) noexcept
{
   const double am = dabs(m);
   const double mtol = tol * am;

   if (mtol == 0.0 && am != 0.0 && tol != 0.0)
      return am <= tol;

   return am <= mtol;
}

constexpr double sign(double x) noexcept
{
   return x >= 0.0 ? 1.0 : -1.0;
}

} // anonymous namespace

/**
 * Returns Re(B0(p^2, m1^2, m2^2, q^2)), from hep-ph/9606211
 */
double b0_softsusy(double p2, double m1, double m2, double q2) noexcept
{
   if (m1 > m2) {
      std::swap(m1, m2);
   }

   // protect against infrared divergence
   if (is_zero(p2, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL)) {
      return 0.0;
   }

   // p is not 0
   if (p2 > 1e-10 * m2) {
      const double s = p2 - m2 + m1;
      const std::complex<double> imin(m1, -EPSTOL);
      const std::complex<double> x = std::sqrt(s*s - 4.0 * p2 * imin);
      const std::complex<double> xp  = (s + sign(s)*x) / (2*p2);
      const std::complex<double> xm = imin / (xp*p2);

      return -std::log(p2/q2) - fB(xp, xm);
   }

   if (is_close(m1, m2, EPSTOL)) {
      return -std::log(m1/q2);
   }

   if (m1 < 1.0e-15) {
      return 1.0 - std::log(m2/q2);
   }

   return 1.0 - std::log(m2/q2)
        + m1 * std::log(m2/m1) / (m1 - m2);
}
