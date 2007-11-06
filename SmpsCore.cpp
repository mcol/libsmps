/*
 *  SmpsCore.cpp
 *
 *  Structure for the core data of the problem.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <string.h>
#include <fstream>
#include "smps.h"
#include "Utils.h"


/** Constructor */
SmpsCore::SmpsCore(const char *coreFileName) :
  nRows(0),
  nStages(0) {
  strcpy(coreFile, coreFileName);
}

/** Constructor */
SmpsCore::SmpsCore(const char *coreFileName, const char *timeFileName) :
  nRows(0),
  nStages(0) {
  strcpy(coreFile, coreFileName);
  strcpy(timeFile, timeFileName);
}

/** Count the number of rows declared in the core file */
int SmpsCore::countRows() {

  ifstream core;
  char buffer[SMPS_LINE_MAX];
  char type[SMPS_FIELD_SIZE], name[SMPS_FIELD_SIZE];
  bool foundRows = false;
  int nValuesRead, rv;

  // open the input file
  core.open(coreFile, ifstream::in);
  if (core.fail()) {
    fprintf(stderr, "Error: Could not open file %s.\n", coreFile);
    return ERROR_FILE_NOT_FOUND;
  }

  // read the file
  while (!core.eof()) {

    // read a line from the file
    rv = readSmpsLine(core, buffer);
    if (rv)
      continue;

    nValuesRead = sscanf(buffer, "%s %s\n", type, name);

    if (nValuesRead == 1) {

      if (strcmp(type, "ROWS") == 0)
	foundRows = true;

      if (strcmp(type, "COLUMNS") == 0)
	break;
    }
    else if (nValuesRead == 2) {
      if (foundRows)
	nRows++;
    }

    else {

    }
  }

  // close the input file
  core.close();

  return 0;
}

/** Count the number of stages declared in the time file */
int SmpsCore::countStages() {

  ifstream time;
  char buffer[SMPS_LINE_MAX];
  char col[SMPS_FIELD_SIZE], row[SMPS_FIELD_SIZE], per[SMPS_FIELD_SIZE];
  int nValuesRead, rv;

  // open the input file
  time.open(timeFile, ifstream::in);
  if (time.fail()) {
    fprintf(stderr, "Error: Could not open file %s.\n", timeFile);
    return ERROR_FILE_NOT_FOUND;
  }

  // read the file
  while (!time.eof()) {

    // read a line from the file
    rv = readSmpsLine(time, buffer);
    if (rv)
      continue;

    nValuesRead = sscanf(buffer, "%s %s %s\n", col, row, per);

    if (nValuesRead == 2) {
      if (strcmp(col, "TIME") != 0) {
	fprintf(stderr, "Something wrong with this line?\n%s\n", buffer);
	return ERROR_TIME_FORMAT;
      }
    }

    else if (nValuesRead == 3) {
      nStages++;
    }
    
    else if (nValuesRead == 1) {

      if (strcmp(col, "ENDATA") == 0) {
	break;
      }
    }
  }

  // close the input file
  time.close();

  return 0;
}
