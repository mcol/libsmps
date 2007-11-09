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
#include "Smps.h"


static int
testCountRows(const char CoreFile[], const int expRows) {

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

static int
testReadCoreFile(const char CoreFile[], const int expStages) {

  int rv;

  char FileName[SMPS_FILENAME_MAX];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, CoreFile);

  SmpsCore core;
  rv = core.readCoreFile(FileName);
  rv = checkEqual(rv, 0, "rv", CoreFile);

  return rv;
}

static int
testReadTimeFile(const char TimeFile[], const int expStages) {

  int rv;

  char FileName[SMPS_FILENAME_MAX];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, TimeFile);

  SmpsCore core;
  rv = core.readTimeFile(FileName);
  if (rv) {
    return TEST_ERROR;
  }

  rv = checkEqual(core.getPeriods(), expStages, "nPeriods", TimeFile);

  return rv;
}

/**
 *  Call the test function for different inputs.
 */
int unitSmpsCore(void) {

  setFamily("testCountRows");
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

  setFamily("testReadCoreFile");
  testReadCoreFile("trivial.cor", 4);
  testReadCoreFile("LandS.cor", 10);
  testReadCoreFile("gbd.cor", 10);
  testReadCoreFile("smallnet.cor", 12);
  testReadCoreFile("mod2-2.cor", 17);
  testReadCoreFile("ssn.cor", 177);
  testReadCoreFile("stocfor1.cor", 118);
  testReadCoreFile("minoux.cor", 118);
  testReadCoreFile("sgpf5y3.cor", 189);
  testReadCoreFile("fxm.cor", 331);
  testReadCoreFile("jll_gva.cor", 350);
  testReadCoreFile("pltexpA2.cor", 167);
  testReadCoreFile("pltexpA7.cor", 687);
  testReadCoreFile("stormG2.cor", 714);
  testReadCoreFile("T1mgnB.cor", 2558);

  setFamily("testReadTimeFile");
  testReadTimeFile("trivial.tim", 2);
  testReadTimeFile("LandS.tim", 2);
  testReadTimeFile("gbd.tim", 2);
  testReadTimeFile("smallnet.tim", 2);
  testReadTimeFile("mod2-2.tim", 2);
  testReadTimeFile("sslp_10_50_100.tim", 2);
  testReadTimeFile("ssn.tim", 2);
  testReadTimeFile("stocfor1.tim", 2);
  testReadTimeFile("minoux.tim", 2);
  testReadTimeFile("sgpf5y3.tim", 3);
  testReadTimeFile("fxm2.tim", 2);
  testReadTimeFile("fxm4.tim", 4);
  testReadTimeFile("jll_gva.tim", 2);
  testReadTimeFile("pltexpA2.tim", 2);
  testReadTimeFile("pltexpA4.tim", 4);
  testReadTimeFile("pltexpA6.tim", 6);
  testReadTimeFile("pltexpA7.tim", 7);
  testReadTimeFile("stormG2.tim", 2);
  testReadTimeFile("T1mgnB.tim", 2);
  testReadTimeFile("timeWAT_10.I", 10);

  return 0;
}
