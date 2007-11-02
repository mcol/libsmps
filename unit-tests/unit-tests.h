/*
 *  unit-tests.h
 *
 *  Common declarations for the unit tests.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#ifndef UNIT_TESTS_H
#define UNIT_TESTS_H

#include <iostream>
using namespace std;

/** Path to the smps test problems directory */
#define SMPS_PATH "/home/mcolombo/oops/smps/testproblems/"

/** Number of test executed */
extern int nTests;

/** Number of test failed */
extern int nFails;

/** Outcomes of the unit tests */
enum TestResults {
  TEST_ERROR = -1,
  TEST_SUCCESS = 0,
  TEST_FAILURE,
  LAST_RESULT
};

/** Call a unit test suite and print a summary of results. */
int callSuite(int (*testSuite)(void), const char name[]);

/** Check equality between two values. */
template<class T>
int checkEqual(const T& value, const T& expValue,
	       const char valueName[], const char testName[] = NULL) {

  nTests++;
  if (value != expValue) {
    if (testName)
      printf("* Testing: %s\n", testName);
    cout << " | FAIL: " << valueName << ": " << value
	 << " (exp: " << expValue << ")\n" << endl;
    nFails++;
    return TEST_FAILURE;
  }

  return TEST_SUCCESS;
}

extern int
unitGetScenarioLength(void);

extern int
unitCountRows(void);

extern int
unitGetStages(void);

extern int
unitTokenizer(void);

#endif /* UNIT_TESTS_H */
