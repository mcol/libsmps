/*
 *  test-SmpsCore.cpp
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
#include "smps.h"

static int
testCountRows(const char CoreFile[], const int expRows);

static int
testGetStages(const char TimeFile[], const int expStages);

/**
 *  Call the test function for different inputs.
 */
int unitCountRows(void) {

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

  return 0;
}

/**
 *  Call the test function for different inputs.
 */
int unitGetStages(void) {

  testGetStages("trivial.tim", 2);
  testGetStages("LandS.tim", 2);
  testGetStages("gbd.tim", 2);
  testGetStages("smallnet.tim", 2);
  testGetStages("mod2-2.tim", 2);
  testGetStages("sslp_10_50_100.tim", 2);
  testGetStages("ssn.tim", 2);
  testGetStages("stocfor1.tim", 2);
  testGetStages("minoux.tim", 2);
  testGetStages("sgpf5y3.tim", 3);
  testGetStages("fxm2.tim", 2);
  testGetStages("fxm4.tim", 4);
  testGetStages("jll_gva.tim", 2);
  testGetStages("pltexpA2.tim", 2);
  testGetStages("pltexpA4.tim", 4);
  testGetStages("pltexpA7.tim", 7);
  testGetStages("stormG2.tim", 2);
  testGetStages("T1mgnB.tim", 2);

  return 0;
}

int testCountRows(const char CoreFile[], const int expRows) {

  int rv;

  char FileName[SMPS_FILENAME_MAX];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, CoreFile);

  /* scan the stochastic file for the number of scenarios and the length */
  SmpsCore core(FileName);
  rv = core.countRows();
  if (rv) {
    return TEST_ERROR;
  }

  rv = checkEqual(core.getRows(), expRows, "nRows", CoreFile);

  return rv;
}

int testGetStages(const char TimeFile[], const int expStages) {

  int rv;

  char FileName[SMPS_FILENAME_MAX];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, TimeFile);

  /* scan the stochastic file for the number of scenarios and the length */
  SmpsCore core("", FileName);
  rv = core.countStages();
  if (rv) {
    return TEST_ERROR;
  }

  rv = checkEqual(core.getStages(), expStages, "nStages", TimeFile);

  return rv;
}
