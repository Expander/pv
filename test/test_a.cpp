#include "pv.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <string>


TEST_CASE("A0")
{
   const auto a0 = [] (double x, double q) {
      return pv::rea0(x, q);
   };

   const std::string filename =
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "A0.txt";
   const auto data = read_from_file(filename);

   for (const auto& v: data) {
      if (v.size() < 3) {
         continue;
      }

      const double x = v.at(0);
      const double q = v.at(1);
      const double a = v.at(2);

      INFO("a0(x = " << x << ", q = " << q << ") = " << a);
      CHECK_CLOSE(a0(x, q), a, 1e-15);
   }
}
