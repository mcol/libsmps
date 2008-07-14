/*
 *  main.cpp
 *
 *  Driver for the SMPS interface to Glpk.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdlib.h>
#include "SmpsGlpk.h"


/** Driver routine for the SMPS interface to Glpk */
int main(const int argc, const char *argv[]) {

  // parse the command line options
  OptionsGlpk opt(argc, argv);
  int rv = opt.parse();
  if (rv)
    return 1;

  // create an object for the problem
  SmpsGlpk data(opt.smpsFile());

  // read the smps files
  rv = data.read();
  if (rv)
    goto TERMINATE;

  // solve the problem
  rv = data.solve(opt);

 TERMINATE:

  return rv;
}
