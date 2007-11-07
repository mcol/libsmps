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

int nTests;
int nFails;

/** Main driver. */
int main(void) {

  callSuite(unitReadSmpsFile, "readSmpsFile");
  callSuite(unitCountRows, "countRows");
  callSuite(unitReadTimeFile, "readTimeFile");
  callSuite(unitReadStocFile, "readStocFile");
  callSuite(unitTokenizer, "Tokenizer");
  callSuite(unitEncodeName, "encodeName");

  return 0;
}

/** Call a unit test suite and print a summary of results. */
int callSuite(int (*testSuite)(void), const char name[]) {

  nTests = 0;
  nFails = 0;
  printf("\n\n *** Tests for %s ***\n", name);

  testSuite();

  printf("\nTotal tests: %d\n", nTests);
  printf("Total fails: %d\n", nFails);

  return 0;
}
