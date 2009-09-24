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

/** Print a separator with the name for a family of tests */
void setFamily(const char *testFamily);

/** Check equality between two values. */
template<class T>
int checkEqual(const T& value, const T& expValue,
	       const string& valueName, const string& testName = "") {

  nTests++;
  if (value != expValue) {
    if (testName != "")
      cout << " | Testing: " << testName << "\n";
    cout << " | FAIL: " << valueName << ": " << value
	 << " (exp: " << expValue << ")\n\n";
    nFails++;
    return TEST_FAILURE;
  }

  return TEST_SUCCESS;
}

#endif /* UNIT_TESTS_H */
