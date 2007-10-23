/*
 *  test-SmpsCore.c
 *
 *  Unit tests for the SmpsCore class.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdio.h>
#include <string.h>
#include "unit-tests.h"
#include "../smps.h"

static int
testCountRows(const char StocFile[], const int expm);

static int nTests = 0;
static int nFails = 0;

/**
 *  Call the testCountRows function and print a summary of results.
 */
int unitCountRows(void) {

  int failures = 0;

  testCountRows("trivial.cor", 4);
  testCountRows("LandS.cor", 10);
  testCountRows("gbd.cor", 10);
  testCountRows("smallnet.cor", 12);
  testCountRows("mod2-2.cor", 17);
  testCountRows("sslp_10_50_100.cor", 62);
  testCountRows("ssn.cor", 177);
  testCountRows("stocfor1.cor", 118);
  testCountRows("minoux.cor", 118);
  testCountRows("sgpf5y3.cor", 189);
  testCountRows("fxm.cor", 331);
  testCountRows("jll_gva.cor", 350);
  testCountRows("pltexpA2.cor", 167);
  testCountRows("pltexpA7.cor", 687);
  testCountRows("stormG2.cor", 714);
  testCountRows("T1mgnB.cor", 2558);

  printf("\nTotal tests: %d\n", nTests);
  printf("Total fails: %d\n", nFails);

  if (nFails == 0 && nTests > 0)
    printf("\n *** ALL TESTS PASS! ***\n");

  return failures;
}

int
testCountRows(const char CoreFile[], const int expRows) {

  int rv;

  char FileName[100];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, CoreFile);

  printf("* Testing: %s\n", CoreFile);

  /* scan the stochastic file for the number of scenarios and the length */
  SmpsCore core(FileName);
  rv = core.countRows();
  if (rv)
    return rv;

  int nRows = core.getRows();
  nTests++;
  if (nRows != expRows) {
    printf(" | FAIL: m: %d (exp: %d)\n", nRows, expRows);
    nFails++;
  }

  return 0;
}
