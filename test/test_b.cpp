#include "pv.hpp"
#include "read_data.hpp"
#include <cmath>
#include <cstdio>

unsigned failed = 0;
unsigned passed = 0;


#define CHECK_EQUAL(a, b)                                       \
   do {                                                         \
      if ((a) == (b)) {                                         \
         passed++;                                              \
      } else {                                                  \
         std::printf("Error: expressions are not equal in line %d: %.17g != %.17g\n", __LINE__, (a), (b));     \
         failed++;                                              \
      }                                                         \
   } while (false);


#define CHECK_CLOSE(a, b, eps)                                  \
   do {                                                         \
      const double diff = std::abs((a) - (b));                  \
      const double mmax = std::max(std::abs(a), std::abs(b));   \
      if (diff <= eps*mmax) {                                   \
         passed++;                                              \
      } else {                                                  \
         std::printf("Error: expressions are not close in line %d: %.17g != %.17g (relative deviation: %g)\n", __LINE__, (a), (b), diff/mmax); \
         failed++;                                              \
      }                                                         \
   } while (false);


template <typename Fn>
void test_b0(Fn fn, double eps)
{
   const std::string filename =
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "B0.txt";
   const auto data = read_from_file(filename);

   for (const auto& v: data) {
      if (v.size() < 6) {
         continue;
      }

      const double p2 = v.at(0);
      const double m1 = v.at(1);
      const double m2 = v.at(2);
      const double q2 = v.at(3);
      const double re = v.at(4);
      // const double im = v.at(5);

      CHECK_CLOSE(fn(p2, m1, m2, q2), re, eps);
   }

   const double p2 = 91*91;
   const double q2 = 100*100;

   // test special cases
   CHECK_CLOSE(fn(2, 3, 4, 1), -1.14772143349266490, eps);
   CHECK_CLOSE(fn(1, 1, 2, 1), -0.263943507355163, eps);
   // CHECK_EQUAL(fn(0, 0, 0, q2), 0.0); // IR divergent case
   CHECK_CLOSE(fn(p2, 0, 0, q2), 2 - log(p2/q2), eps);
   CHECK_CLOSE(fn(0, p2, 0, q2), 1 - log(p2/q2), eps);
   CHECK_CLOSE(fn(p2, p2, 0, q2), 2 - log(p2/q2), eps);
   CHECK_CLOSE(fn(0, p2, p2, q2), -std::log(p2/q2), eps);

   // symmetries
   CHECK_EQUAL(fn(0, 0, p2, q2), fn(0, p2, 0, q2));
   CHECK_EQUAL(fn(p2, 0, p2, q2), fn(p2, p2, 0, q2));
}


template <typename Fn>
void test_b0_limits(Fn fn, double eps)
{
   const std::string filename =
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "B0_limits.txt";
   const auto data = read_from_file(filename);

   for (const auto& v: data) {
      if (v.size() < 5) {
         continue;
      }

      const double p2 = v.at(0);
      const double m1 = v.at(1);
      const double m2 = v.at(2);
      const double q2 = v.at(3);
      const double re = v.at(4);

      CHECK_CLOSE(fn(p2, m1, m2, q2), re, eps);
   }
}


template <typename Fn>
void test_b0_xx(Fn fn, double eps)
{
   const std::string filename =
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "B0_xx.txt";
   const auto data = read_from_file(filename);

   for (const auto& v: data) {
      if (v.size() < 5) {
         continue;
      }

      const double p2 = v.at(0);
      const double m1 = v.at(1);
      const double m2 = v.at(2);
      const double q2 = v.at(3);
      const double re = v.at(4);

      CHECK_CLOSE(fn(p2, m1, m2, q2), re, eps);
   }
}


void test_b0()
{
   const auto b0 = [] (double s, double x, double y, double q) {
      return pv::b0(s, x, y, q);
   };
   const auto degrassi = [] (double s, double x, double y, double q) {
      return b0_degrassi(s, x, y, q);
   };
   const auto softsusy = [] (double s, double x, double y, double q) {
      return pv::b0_softsusy(s, x, y, q);
   };

   test_b0(b0, 3e-11);
   test_b0(degrassi, 3e-11);
   test_b0(softsusy, 1e-05);
   // test_b0_limits(degrassi, 1e-02);
   // test_b0_limits(softsusy, 1e-02);
   // test_b0_xx(degrassi, 1e-02);
   // test_b0_xx(softsusy, 1e-02);
}


int main()
{
   test_b0();

   std::printf("Passed tests: [%u/%u]\n", passed, passed + failed);

   return failed != 0;
}
