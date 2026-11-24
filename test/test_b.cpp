#include "alt.hpp"
#include "pv.hpp"
#include "read_data.hpp"
#include "test.hpp"


TEST_CASE("B0 grid")
{
   const std::string filename =
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "B0.txt";
   const auto data = read_from_file(filename);

   for (const auto& v: data) {
      if (v.size() < 6) {
         continue;
      }

      const double s = v.at(0);
      const double x = v.at(1);
      const double y = v.at(2);
      const double q = v.at(3);
      const double b = v.at(4);
      // const double im = v.at(5);

      INFO("b0(s = " << s << ", x = " << x << ", y = " << y << ", q = " << q << ") = " << b);
      CHECK_CLOSE(pv::reb0(s, x, y, q), b, 3e-11);
      CHECK_CLOSE(b0_degrassi(s, x, y, q), b, 3e-11);
      CHECK_CLOSE(pv::b0_gpt_astra(s, x, y, q), b, 5e-08);
      CHECK_CLOSE(pv::b0_softsusy(s, x, y, q), b, 1e-05);
   }
}


TEST_CASE("B0 special cases")
{
   const double s = 91*91;
   const double q = 100*100;

   // test special cases
   CHECK_CLOSE(pv::reb0(2, 3, 4, 1), -1.14772143349266490, 3e-11);
   CHECK_CLOSE(pv::reb0(1, 1, 2, 1), -0.263943507355163, 3e-11);
   CHECK_EQUAL(pv::reb0(0, 0, 0, q), 0.0); // IR divergent case
   CHECK_CLOSE(pv::reb0(s, 0, 0, q), 2 - log(s/q), 3e-11);
   CHECK_CLOSE(pv::reb0(0, s, 0, q), 1 - log(s/q), 3e-11);
   CHECK_CLOSE(pv::reb0(s, s, 0, q), 2 - log(s/q), 3e-11);
   CHECK_CLOSE(pv::reb0(0, s, s, q), -std::log(s/q), 3e-11);

   // symmetries
   CHECK_EQUAL(pv::reb0(0, 0, s, q), pv::reb0(0, s, 0, q));
   CHECK_EQUAL(pv::reb0(s, 0, s, q), pv::reb0(s, s, 0, q));

   // renormalization scale <= 0
   if (std::numeric_limits<double>::has_infinity) {
      CHECK(std::isinf(pv::reb0(0, 0, 0, 0)));
      CHECK(std::isinf(pv::reb0(1, 2, 3, 0)));
      CHECK(std::isinf(pv::reb0(1, 2, 3, -1)));
   }

   // IR divergence
   CHECK(pv::reb0(0, 0, 0, 1) == 0);

   // small input
   {
      constexpr double min = std::numeric_limits<double>::min();
      constexpr double pi = 3.1415926535897932;
      constexpr double q = 1;
      CHECK(pv::reb0(min, min, min, q) == 2 - pi/std::sqrt(3) - std::log(min/q));
   }
}


TEST_CASE("B0 limits")
{
   const std::string filename =
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "B0_limits.txt";
   const auto data = read_from_file(filename);

   for (const auto& v: data) {
      if (v.size() < 5) {
         continue;
      }

      const double s = v.at(0);
      const double x = v.at(1);
      const double y = v.at(2);
      const double q = v.at(3);
      const double b = v.at(4);

      INFO("b0(s = " << s << ", x = " << x << ", y = " << y << ", q = " << q << ") = " << b);
      CHECK_CLOSE(pv::reb0(s, x, y, q), b, 1e-08);
      CHECK_CLOSE(pv::b0_gpt_astra(s, x, y, q), b, 1e-12);
   }
}


TEST_CASE("B0 equal masses")
{
   const std::string filename =
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "B0_xx.txt";
   const auto data = read_from_file(filename);

   for (const auto& v: data) {
      if (v.size() < 5) {
         continue;
      }

      const double s = v.at(0);
      const double x = v.at(1);
      const double y = v.at(2);
      const double q = v.at(3);
      const double b = v.at(4);

      INFO("b0(s = " << s << ", x = " << x << ", y = " << y << ", q = " << q << ") = " << b);
      CHECK_CLOSE(pv::reb0(s, x, y, q), b, 1e-12);
   }
}
