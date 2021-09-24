#include "pv.hpp"
#include "read_data.hpp"
#include <cstdio>

unsigned failed = 0;
unsigned passed = 0;


#define CHECK_EQ(a, b)                                          \
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


void test_b0()
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

      CHECK_CLOSE(b0_degrassi(p2, m1, m2, q2), re, 2e-11);
      CHECK_CLOSE(b0_softsusy(p2, m1, m2, q2), re, 1e-05);
   }
}


int main()
{
   test_b0();

   std::printf("Passed tests: [%u/%u]\n", passed, passed + failed);

   return failed != 0;
}
