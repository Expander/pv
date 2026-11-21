#pragma once

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#define CHECK_CLOSE(a, b, eps) CHECK((a) == doctest::Approx(b).epsilon(eps));

#define CHECK_EQUAL(a, b) CHECK((a) == (b))
