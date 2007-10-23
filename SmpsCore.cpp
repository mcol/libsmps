#include <stdio.h>
#include <string.h>
#include "smps.h"


/** Constructor */
SmpsCore::SmpsCore(const char *coreFileName) {
  
  strcpy(coreFile, coreFileName);
  nRows = 0;
}

/** Retrieve the number of rows in the core file */
int SmpsCore::getRows() {
  return nRows;
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
