/*
 *  unit-tests.cpp
 *
 *  Main driver for the unit tests.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include "unit-tests.h"


extern int
unitSmps(void);

extern int
unitSmpsCore(void);

extern int
unitSmpsTree(void);

extern int
unitTokenizer(void);

extern int
unitUtils(void);

int nTests;
int nFails;

/** Main driver. */
int main(void) {

  callSuite(unitSmps,      "Smps");
  callSuite(unitSmpsCore,  "SmpsCore");
  callSuite(unitSmpsTree,  "SmpsTree");
  callSuite(unitTokenizer, "Tokenizer");
  callSuite(unitUtils,     "Utils");

  return 0;
}

/** Call a unit test suite and print a summary of results. */
int callSuite(int (*testSuite)(void), const char name[]) {

  nTests = 0;
  nFails = 0;
  printf("\n *** Tests for %s ***\n", name);

  testSuite();

  printf(" === Total tests: %d\n", nTests);
  printf(" === Total fails: %d\n", nFails);

  return 0;
}

/** Print a separator with the name for a family of tests */
void setFamily(const char *testFamily) {
  printf("   = %s\n", testFamily);
}
