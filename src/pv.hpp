#pragma once

extern "C" {

// real part of B0, from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104
double b0_degrassi(double p2, double m1, double m2, double q2);

}

namespace pv {

// optimized implementation from [hep-ph/9606211]
double b0_softsusy(double p2, double m1, double m2, double q2) noexcept;

// ultimate implementation
double b0(double s, double x, double y, double q) noexcept;

} // namespace pv
