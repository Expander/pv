#pragma once


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
