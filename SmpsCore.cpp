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
#include "smps.h"


/** Constructor */
SmpsCore::SmpsCore(const char *coreFileName) {
  
  strcpy(coreFile, coreFileName);
  nRows = 0;
  nStages = 0;
}

/** Constructor */
SmpsCore::SmpsCore(const char *coreFileName, const char *timeFileName) {
  
  strcpy(coreFile, coreFileName);
  strcpy(timeFile, timeFileName);
  nRows = 0;
  nStages = 0;
}

/** Retrieve the number of rows in the core file */
int SmpsCore::getRows() const {
  return nRows;
}

/** Retrieve the number of stages in the problem */
int SmpsCore::getStages() const {
  return nStages;
}

/** Count the number of rows declared in the core file */
int SmpsCore::countRows() {

  char buffer[LINE_MAX];
  char type[SMPS_FIELD_SIZE], name[SMPS_FIELD_SIZE];
  bool foundRows = false;

  int nValuesRead;

  // open the input file
  FILE *core = fopen(coreFile, "r");
  if (!core) {
    fprintf(stderr, "Error: Could not open file %s.\n", coreFile);
    return ERROR_FILE_NOT_FOUND;
  }

  // read all lines of the file
  while (fgets(buffer, LINE_MAX, core) != NULL) {

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
  fclose(core);

  return 0;
}

/** Count the number of stages declared in the time file */
int SmpsCore::countStages() {

  char buffer[LINE_MAX];
  char col[SMPS_FIELD_SIZE], row[SMPS_FIELD_SIZE], per[SMPS_FIELD_SIZE];
  int nValuesRead;

  // open the input file
  FILE *time = fopen(timeFile, "r");
  if (!time) {
    fprintf(stderr, "Error: Could not open file %s.\n", timeFile);
    return ERROR_FILE_NOT_FOUND;
  }

  // read all lines of the file
  while (fgets(buffer, LINE_MAX, time) != NULL) {

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
  fclose(time);

  return 0;
}
