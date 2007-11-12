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

  Smps tree(FileName);
  tree.readSmpsFile();
  tree.readCoreFile();
  tree.readTimeFile();
  rv = tree.readStocFile();
  if (rv) {
    return TEST_ERROR;
  }

  rv = checkEqual(tree.getMaxScens(), expScens, "maxScens", StocFile);
  rv = checkEqual(tree.getMaxNodes(), expNodes, "maxNodes", StocFile);
  rv = checkEqual(tree.getMaxReals(), expReals, "maxReals", StocFile);

  return rv;
}

/**
 *  Call the test function for different inputs.
 */
int unitSmpsTree(void) {

#ifndef SKIP_TEST_INDEP
  setFamily("testReadStocFile Indep");
  testReadStocFile("trivial.smps", 4, 5, 9);
  testReadStocFile("trivial2.smps", 2, 3, 3);
  testReadStocFile("trivial4.smps", 4, 5, 5);
  testReadStocFile("trivial5.smps", 2, 3, 3);
  testReadStocFile("trivial-mtx.smps", 2, 3, 5);
  testReadStocFile("fxm2_1.smps", 1, 2, 2);
  testReadStocFile("fxm2_2.smps", 2, 3, 3);
  testReadStocFile("fxm2_6.smps", 6, 7, 7);
  testReadStocFile("fxm3_1.smps", 1, 3, 3);
  testReadStocFile("fxm3_2.smps", 4, 7, 9);
  testReadStocFile("fxm3_6.smps", 36, 43, 73);
  testReadStocFile("fxm3_16.smps", 256, 337, 769);
  testReadStocFile("fxm4_6.smps", 216, 259, 216);
#endif /* SKIP_TEST_INDEP */

#ifndef SKIP_TEST_BLOCK
  setFamily("testReadStocFile Blocks");
  testReadStocFile("trivial-blk.smps", 2, 3, 5);
  testReadStocFile("LandS-blocks.smps", 3, 4, 4);
  testReadStocFile("mod2-2.smps", 10, 11, 101);
  testReadStocFile("stocfor1.smps", 1, 2, 15);
  testReadStocFile("stocfor2.smps", 64, 65, 5277);
  testReadStocFile("ft3-small.smps", 2, 3, 65);
  testReadStocFile("pltexpA2_6.smps", 6, 7, 43);
  testReadStocFile("pltexpA3_6.smps", 36, 43, 505);
  testReadStocFile("pltexpA4_6.smps", 216, 259, 4537);
  testReadStocFile("pltexpA6_6.smps", 7776, 9331, 272161);
  testReadStocFile("storm4.smps", 4, 5, 311);
  testReadStocFile("storm8.smps", 8, 9, 937);
  testReadStocFile("minoux-100-0.5.smps", 100, 101, 6601);
#endif /* SKIP_TEST_BLOCK */

  return 0;
}
