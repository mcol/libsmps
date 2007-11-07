/*
 *  test-SmpsTree.cpp
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
#include "Smps.h"


static int
testReadStocFile(const char StocFile[], const int expScens,
		 const int expNodes, const int expReals) {

  int rv;

  char FileName[SMPS_FILENAME_MAX];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, StocFile);

  SmpsTree tree(FileName);
  rv = tree.readStocFile(FileName);
  if (rv) {
    return TEST_ERROR;
  }

  rv = checkEqual(tree.getMaxScens(), expScens, "maxScens", StocFile);
  rv = checkEqual(tree.getMaxNodes(), expNodes, "maxNodes", StocFile);
  //  rv = checkEqual(tree.getMaxReals(), expReals, "maxReals", StocFile);

  return rv;
}

/**
 *  Call the test function for different inputs.
 */
int unitReadStocFile(void) {

#ifndef SKIP_TEST_INDEP
  testReadStocFile("trivial.sto", 4, 5, 4);
  testReadStocFile("trivial2.sto", 2, 3, 2);
  testReadStocFile("trivial4.sto", 4, 5, 4);
  testReadStocFile("trivial5.sto", 2, 3, 2);
  testReadStocFile("trivial-mtx.sto", 2, 3, 3);
  testReadStocFile("fxm2_1.sto", 1, 2, 1);
  testReadStocFile("fxm2_2.sto", 2, 3, 2);
  testReadStocFile("fxm2_6.sto", 6, 7, 6);
  testReadStocFile("fxm3_1.sto", 1, 3, 2);
  testReadStocFile("fxm3_2.sto", 4, 7, 4);
  testReadStocFile("fxm3_6.sto", 36, 43, 12);
  testReadStocFile("fxm4_6.sto", 216, 259, 18);
#endif /* SKIP_TEST_INDEP */

#ifndef SKIP_TEST_BLOCK
  testReadStocFile("trivial-blk.sto", 2, 3, 0);
  testReadStocFile("LandS_blocks.sto", 3, 4, 0);
  testReadStocFile("mod2-2.sto", 10, 11, 0);
  testReadStocFile("stocfor1.sto", 1, 2, 0);
  testReadStocFile("stocfor2.sto", 64, 65, 0);
  testReadStocFile("ft3-small.sto", 2, 3, 0);
  testReadStocFile("pltexpA2_6.sto", 6, 7, 0);
  testReadStocFile("pltexpA3_6.sto", 36, 43, 0);
  testReadStocFile("pltexpA4_6.sto", 216, 259, 0);
  testReadStocFile("pltexpA6_6.sto", 7776, 9331, 0);
  testReadStocFile("stormG2_4.sto", 4, 5, 0);
  testReadStocFile("stormG2_8.sto", 8, 9, 0);
  testReadStocFile("minoux-100-0.5.sto", 100, 101, 0);
#endif /* SKIP_TEST_BLOCK */

  return 0;
}
