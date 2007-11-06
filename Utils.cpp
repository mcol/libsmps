/*
 *  Utils.cpp
 *
 *  Various utility functions.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <string.h>
#include "Utils.h"


/** Encode the name into an integer */
int encodeName(const char* name, const int nRows) {

  int code = 0;
  int nameLength = strlen(name);

  // compute a numerical encoding of the given name
  for (int i = 0; i < nameLength; ++i)
    code += name[i] * (i + 1);

  // restrict the value of the code in [1, nRows]
  code = (code % nRows) + 1;

  return code;
}

/** Determine the index such that names[index] == name */
int lookupCode(const char **names, const char *name, int nRows,
	       int *header, int *links) {

  int kCode, index;

  // get the code for name
  kCode = encodeName(name, nRows);

  // look up the initial index of the code
  index = header[kCode];

  // find the index of the code by walking through the linked list of codes
  for (int i = 0; i < nRows; i++) {

    if (index == 0)
      break;

    if (strcmp(names[index], name) == 0)
      break;

    index = links[index];
  }

  return index;
}


/** Read a line from an Smps file */
int readSmpsLine(ifstream &file, char *buffer) {

  // read a line from the file
  file.getline(buffer, SMPS_LINE_MAX);

  // skip the asterisk lines
  if (buffer[0] == '*')
    return 1;

#ifdef DEBUG_SMPS_BUFFER
  printf(buffer);
#endif

  return 0;
}
