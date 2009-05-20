/*
 *  main.cpp
 *
 *  Driver for the SMPS interface to Oops.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include "SmpsOops.h"


FILE *printout = stdout;
FILE *globlog;


/** Driver routine for the SMPS interface to OOPS */
int main(const int argc, const char *argv[]) {

  // initialise the parallel stuff
  InitLippPar(argc, (char **) argv);

  // parse the command line options
  OptionsOops opt(argc, argv);
  int rv = opt.parse();
  if (rv)
    return 1;

  // create an object for the problem
  SmpsOops data(opt.smpsFile(), opt.cutoffLevel());

  // read the smps files
  rv = data.read();
  if (rv)
    return 1;

  HopdmOptions hopdmOpts = HopdmOptions();

  // warmstart case
  if (opt.useWarmstart()) {

    // select the warmstart strategy
    if (opt.useReduction())
      rv = data.reduceTree(opt.useReduction());
    else if (opt.useAggregation())
      rv = data.aggregateStages(opt.useAggregation());

    if (rv)
      goto LEAVE_WARMSTART;

    // solve the reduced problem
    rv = data.solveReduced(opt, hopdmOpts);
    if (rv)
      goto TERMINATE;
  }

 LEAVE_WARMSTART:

  // solve the problem
  rv = data.solve(opt, hopdmOpts);

 TERMINATE:

#ifdef WITH_MPI
  MPI_Finalize();
#endif

  return rv;
}
