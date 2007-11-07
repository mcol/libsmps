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
testReadSmpsFile(const char SmpsFile[]) {

  int rv;

  char FileName[SMPS_FILENAME_MAX];
  char SmpsPath[] = SMPS_PATH;
  strcpy(FileName, SmpsPath);
  strcat(FileName, SmpsFile);

  Smps smps(FileName);
  rv = smps.readSmpsFile(FileName);

  // here we only test that the operation was successful
  rv = checkEqual(rv, 0, "rv", SmpsFile);

  return rv;
}

/**
 *  Call the test function for different inputs.
 */
int unitSmps(void) {

  testReadSmpsFile("fxm3_2.smps");
  testReadSmpsFile("gbd.smps");
  testReadSmpsFile("LandS.smps");
  testReadSmpsFile("mod2-2.smps");
  testReadSmpsFile("sgpf5y3.smps");
  testReadSmpsFile("smallnet.smps");
  testReadSmpsFile("ssn.smps");
  testReadSmpsFile("stocfor1.smps");
  testReadSmpsFile("storm8.smps");
  testReadSmpsFile("trivial.smps");

  return 0;
}
