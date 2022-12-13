#include "alt.hpp"
#include "bench.hpp"
#include "pv.hpp"
#include <array>
#include <string>
#include <iomanip>
#include <iostream>
#include <vector>

template <typename T, typename Fn>
void bench_fn(Fn f, const std::vector<T>& values, const std::string& name,
              const std::string& type)
{
   // warm-up
   for (const auto& v: values) {
      pv::bench::do_not_optimize(f(v[0], v[1], v[2], v[3]));
   }

   const auto total_time = pv::bench::time_in_seconds([&] {
         for (const auto& v: values) {
            pv::bench::do_not_optimize(f(v[0], v[1], v[2], v[3]));
         }
      });

   std::cout << std::setw(24) << std::left << name << "type: " << std::setw(16)
             << std::left << type << "time: " << total_time << "s\n";
}

void print_line(char c)
{
   for (int i = 0; i < 60; ++i) {
      std::cout << c;
   }
   std::cout << '\n';
}

void print_headline(const std::string& text)
{
   print_line('=');
   std::cout << text << '\n';
   print_line('=');
}

template<typename T>
void bench(const T& values_d)
{
   print_headline("B0");

   bench_fn([&](double s, double x, double y, double q) {
               return pv::reb0(s,x,y,q);
            },
            values_d, "B0", "real");

   bench_fn([&](double s, double x, double y, double q) {
               return b0_degrassi(s,x,y,q);
            },
            values_d, "Degrassi", "real");

   bench_fn([&](double s, double x, double y, double q) {
               return pv::b0_softsusy(s,x,y,q);
            },
            values_d, "SOFTSUSY", "real");
}

int main()
{
   using pv::bench::generate_random_scalars;

   const std::size_t N = 1'000'000;
   const auto min = 0;
   const auto max = 4;

   const auto s = generate_random_scalars<double>(N, min, max);
   const auto x = generate_random_scalars<double>(N, min, max);
   const auto y = generate_random_scalars<double>(N, min, max);
   const auto q = generate_random_scalars<double>(N, min, max);

   std::vector<std::array<double, 4>> values(N);

   for (std::size_t i = 0; i < N; ++i) {
      values.at(i).at(0) = s.at(i);
      values.at(i).at(1) = x.at(i);
      values.at(i).at(2) = y.at(i);
      values.at(i).at(3) = q.at(i);
   }

   bench(values);

   return 0;
}
