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


/** Encode the name into an integer */
int encodeName(const char* name, const int m);

/** Determine the index such that names[index] == name */
int lookupCode(const char **names, const char *name, int m,
	       int *header, int *links);

#endif /* _UTILS_H_ */
