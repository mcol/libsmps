/*
 *  Utils.h
 *
 *  Various utility functions.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#ifndef _UTILS_H_
#define _UTILS_H_

#include <fstream>
#include "Smps.h"
using namespace std;


/** Encode the name into an integer */
int encodeName(const char* name, const int m);

/** Determine the index such that names[index] == name */
int lookupCode(const char **names, const char *name, int m,
	       int *header, int *links);

/** Read a line from an Smps file */
int readSmpsLine(ifstream &file, char *buffer);

#endif /* _UTILS_H_ */
