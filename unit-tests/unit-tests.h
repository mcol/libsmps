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

/** Path to the smps test problems directory */
#define SMPS_PATH "/home/mcolombo/oops/smps/testproblems/"

/** Outcomes of the unit tests */
enum TestResults {
  TEST_ERROR = -1,
  TEST_SUCCESS = 0,
  TEST_FAILURE,
  LAST_RESULT
};

/** Call a unit test suite and print a summary of results. */
int callSuite(int (*testSuite)(void), const char name[]);

extern int
unitGetScenarioLength(void);

extern int
unitCountRows(void);

extern int
unitGetStages(void);

extern int
unitTokenizer(void);

#endif /* UNIT_TESTS_H */
