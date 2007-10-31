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

#include <stdio.h>
#include <cstring>
#include "unit-tests.h"

int nTests;
int nFails;

/** Main driver. */
int main(void) {

  callSuite(unitGetScenarioLength, "getScenarioLength");
  callSuite(unitCountRows, "countRows");
  callSuite(unitGetStages, "getStages");
  callSuite(unitTokenizer, "Tokenizer");

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

/** Check numerical equality between two integer values */
int checkEqual(const int value, const int expValue,
	       const char valueName[], const char testName[]) {

  nTests++;
  if (value != expValue) {
    if (testName)
      printf("* Testing: %s\n", testName);
    printf(" | FAIL: %s: %d (exp: %d)\n", valueName, value, expValue);
    nFails++;
    return TEST_FAILURE;
  }

  return TEST_SUCCESS;
}

/** Check numerical equality between two strings */
int checkEqual(const char value[], const char expValue[],
	       const char valueName[], const char testName[]) {

  nTests++;
  if (strcmp(value, expValue) != 0) {
    if (testName)
      printf("* Testing: %s\n", testName);
    printf(" | FAIL: %s: '%s' (exp: '%s')\n", valueName, value, expValue);
    nFails++;
    return TEST_FAILURE;
  }

  return TEST_SUCCESS;
}
