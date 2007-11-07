/*
 *  test-Smps.cpp
 *
 *  Unit tests for the Smps class.
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
testReadTimeFile(const char TimeFile[], const int expStages) {

  int rv;

  char FileName[SMPS_FILENAME_MAX];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, TimeFile);

  Smps smps;
  rv = smps.readTimeFile(FileName);
  if (rv) {
    return TEST_ERROR;
  }

  rv = checkEqual(smps.getPeriods(), expStages, "nPeriods", TimeFile);

  return rv;
}

/**
 *  Call the test function for different inputs.
 */
int unitReadTimeFile(void) {

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
