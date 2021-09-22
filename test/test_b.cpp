#include "pv.h"
#include <cstdio>

unsigned failed = 0;
unsigned passed = 0;

void check(bool cond)
{
   if (cond) {
      passed++;
   } else {
      std::printf("test failed\n");
      failed++;
   }
}

void test_b0()
{
   check(b0_degrassi(0,1,1,1) == 0.0);
}

int main()
{
   test_b0();

   std::printf("Passed tests: [%u/%u]\n", passed, passed + failed);

   return failed != 0;
}
