#include "pv.hpp"
#include "read_data.hpp"
#include "test.hpp"
#include <string>


TEST_CASE("A0")
{
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
      CHECK_CLOSE(pv::rea0(x, q), a, 1e-15);
      CHECK(pv::rea0(x, q) == -pv::rea0(-x, q));
   }
}


TEST_CASE("A0 special cases")
{
   CHECK(!std::isfinite(pv::rea0(0, 0)));
   CHECK(pv::rea0(0, 1) == 0);
}
