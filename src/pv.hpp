#pragma once

extern "C" {

/* real part of B0, from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104 */
double b0_degrassi(double p2, double m1, double m2, double q2);

}

double b0_softsusy(double p, double m1, double m2, double q) noexcept;
