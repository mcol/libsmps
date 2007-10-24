/*
 *  unit-tests.c
 *
 *  Main driver for the unit tests.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdio.h>
#include "unit-tests.h"

int nTests;
int nFails;

/** Main driver. */
int main(void) {

  callSuite(unitGetScenarioLength, "getScenarioLength");
  callSuite(unitCountRows, "countRows");
  callSuite(unitGetStages, "getStages");

  return 0;
}

/** Call a unit test suite and print a summary of results. */
int callSuite(int (*testSuite)(void), const char name[]) {

  nTests = 0;
  nFails = 0;
  printf("\n\n *** Tests for %s ***\n\n", name);

  testSuite();

  printf("\nTotal tests: %d\n", nTests);
  printf("Total fails: %d\n", nFails);

  return 0;
}
