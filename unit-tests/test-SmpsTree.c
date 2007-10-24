/*
 *  test-SmpsTree.c
 *
 *  Unit tests for the SmpsTree class.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdio.h>
#include <string.h>
#include "unit-tests.h"
#include "smps.h"

#define TEST_INDEP 1
#define TEST_BLOCK 1

extern int nTests;
extern int nFails;

static int
testGetScenarioLength(const char StocFile[], const int expScens,
		      const int expNodes, const int expReals);

/**
 *  Call the test function for different inputs.
 */
int unitGetScenarioLength(void) {

#if TEST_INDEP
  testGetScenarioLength("trivial.sto", 4, 5, 4);
  testGetScenarioLength("trivial2.sto", 2, 3, 2);
  testGetScenarioLength("trivial4.sto", 4, 5, 4);
  testGetScenarioLength("trivial5.sto", 2, 3, 2);
  testGetScenarioLength("trivial-mtx.sto", 2, 3, 3);
  testGetScenarioLength("fxm2_1.sto", 1, 2, 1);
  testGetScenarioLength("fxm2_2.sto", 2, 3, 2);
  testGetScenarioLength("fxm2_6.sto", 6, 7, 6);
  testGetScenarioLength("fxm3_1.sto", 1, 3, 2);
  testGetScenarioLength("fxm3_2.sto", 4, 7, 4);
  testGetScenarioLength("fxm3_6.sto", 36, 43, 12);
  testGetScenarioLength("fxm4_6.sto", 216, 259, 18);
#endif /* TEST_INDEP */

#if TEST_BLOCK
  testGetScenarioLength("trivial-blk.sto", 2, 3, 0);
  testGetScenarioLength("LandS_blocks.sto", 3, 4, 0);
  testGetScenarioLength("mod2-2.sto", 10, 11, 0);
  testGetScenarioLength("stocfor1.sto", 1, 2, 0);
  testGetScenarioLength("stocfor2.sto", 64, 65, 0);
  testGetScenarioLength("ft3-small.sto", 2, 3, 0);
  testGetScenarioLength("pltexpA2_6.sto", 6, 7, 0);
  testGetScenarioLength("pltexpA3_6.sto", 36, 43, 0);
  testGetScenarioLength("pltexpA4_6.sto", 216, 259, 0);
  testGetScenarioLength("pltexpA6_6.sto", 7776, 9331, 0);
  testGetScenarioLength("stormG2_4.sto", 4, 5, 0);
  testGetScenarioLength("stormG2_8.sto", 8, 9, 0);
  testGetScenarioLength("minoux-100-0.5.sto", 100, 101, 0);
#endif /* TEST_BLOCK */

  return 0;
}

int
testGetScenarioLength(const char StocFile[], const int expScens,
		      const int expNodes, const int expReals) {

  int rv;

  char FileName[100];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, StocFile);

  printf("* Testing: %s\n", StocFile);

  /* scan the stochastic file for the number of scenarios and the length */
  SmpsTree tree(FileName);
  rv = tree.getScenarioLength();
  if (rv)
    return rv;

  int maxScens = tree.getMaxScens();
  nTests++;

  int maxNodes = tree.getMaxNodes();
  nTests++;

  //  int maxReals = tree.getMaxReals();
  //  nTests++;

  if (maxScens != expScens) {
    printf(" | FAIL: maxScens: %d (exp: %d)\n", maxScens, expScens);
    nFails++;
  }
  if (maxNodes != expNodes) {
    printf(" | FAIL: maxNodes: %d (exp: %d)\n", maxNodes, expNodes);
    nFails++;
  }
  //  if (maxReals != expReals) {
  //    printf(" | FAIL: maxReals: %d (exp: %d)\n", maxReals, expReals);
  //    nFails++;
  //  }

  return 0;
}
