#include "pv.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <string>


template <typename Fn>
void test_a0(Fn fn, double eps)
{
   const std::string filename =
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "A0.txt";
   const auto data = read_from_file(filename);

   for (const auto& v: data) {
      if (v.size() < 3) {
         continue;
      }

      const double m2 = v.at(0);
      const double q2 = v.at(1);
      const double re = v.at(2);

      CHECK_CLOSE(fn(m2, q2), re, eps);
   }
}


void test_a0()
{
   const auto a0 = [] (double x, double q) {
      return pv::rea0(x, q);
   };

   test_a0(a0, 1e-15);
}


int main()
{
   test_a0();

   std::printf("Failed tests: [%u/%u]\n", failed, passed + failed);
   std::printf("Passed tests: [%u/%u]\n", passed, passed + failed);

   return failed != 0;
}
