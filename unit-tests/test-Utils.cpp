/*
 *  test-Utils.cpp
 *
 *  Unit tests for the Utils functions.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include "unit-tests.h"
#include "Utils.h"

static int
testEncodeName(const char name[], const int m, const int expCode) {

  int rv;

  int code = encodeName(name, m);

  rv = checkEqual(code, expCode, "code", name);

  return rv;
}

/**
 *  Call the test function for different inputs.
 */
int unitEncodeName(void) {

  testEncodeName("trivial", 4, 3);
  testEncodeName("LandS", 10, 6);
  testEncodeName("gbd", 10, 10);
  testEncodeName("smallnet", 12, 4);
  testEncodeName("mod2", 17, 16);
  testEncodeName("sslp", 62, 2);
  testEncodeName("ssn", 177, 145);
  testEncodeName("stocfor1", 118, 21);
  testEncodeName("minoux", 118, 39);
  testEncodeName("sgpf5y3", 189, 146);
  testEncodeName("fxm", 331, 8);
  testEncodeName("jll_gva", 350, 129);
  testEncodeName("pltexpA2", 167, 35);
  testEncodeName("pltexpA7", 687, 500);
  testEncodeName("stormG2", 714, 316);
  testEncodeName("T1mgnB", 2558, 1868);
  testEncodeName("", 200, 1);

  return 0;
}
