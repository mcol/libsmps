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
testReadSmpsFile(const string SmpsFile) {

  int rv;
  const string FileName = SMPS_PATH + SmpsFile;
  Smps smps(FileName);

  rv = smps.readSmpsFile(FileName);
  rv = checkEqual(rv, 0, "readSmpsFile rv", SmpsFile);
  if (rv)
    return rv;

  rv = smps.readCoreFile();
  rv = checkEqual(rv, 0, "readCoreFile rv", SmpsFile);

  rv = smps.readTimeFile();
  rv = checkEqual(rv, 0, "readTimeFile rv", SmpsFile);

  rv = smps.readStocFile();
  rv = checkEqual(rv, 0, "readStocFile rv", SmpsFile);

  return rv;
}

/**
 *  Call the test function for different inputs.
 */
int unitSmps(void) {

  setFamily("testReadSmpsFile");
  testReadSmpsFile("fxm3_2.smps");
  testReadSmpsFile("gbd.smps");
  testReadSmpsFile("LandS.smps");
  testReadSmpsFile("mod2-2.smps");
  //  testReadSmpsFile("sgpf5y3.smps");
  testReadSmpsFile("smallnet.smps");
  //  testReadSmpsFile("ssn.smps");
  testReadSmpsFile("stocfor1.smps");
  testReadSmpsFile("storm8.smps");
  testReadSmpsFile("trivial.smps");

  return 0;
}
