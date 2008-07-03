/*
 *  main.cpp
 *
 *  Driver for the SMPS interface to lpsolve.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdlib.h>
#include "SmpsLpsolve.h"


/** Driver routine for the SMPS interface to lpsolve */
int main(const int argc, const char *argv[]) {

  // parse the command line options
  OptionsLpsolve opt(argc, argv);
  int rv = opt.parse();
  if (rv)
    return 1;

  // create an object for the problem
  SmpsLpsolve data(opt.smpsFile());

  // read the smps files
  rv = data.read();
  if (rv)
    goto TERMINATE;

  // solve the problem
  rv = data.solve(opt);

 TERMINATE:

  return rv;
}
