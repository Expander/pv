#include "pv.h"
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


void test_b0()
{
   CHECK_EQ(b0_degrassi(0,0,1,1), 1.0);
   CHECK_EQ(b0_degrassi(0,1,1,1), 0.0);
}


int main()
{
   test_b0();

   std::printf("Passed tests: [%u/%u]\n", passed, passed + failed);

   return failed != 0;
}
