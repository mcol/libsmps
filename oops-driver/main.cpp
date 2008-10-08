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


FILE *printout;
static void setupOutputFile(OptionsOops &opt);


/** Driver routine for the SMPS interface to OOPS */
int main(const int argc, const char *argv[]) {

  // initialise the parallel stuff
  InitLippPar(argc, (char **) argv);

  // parse the command line options
  OptionsOops opt(argc, argv);
  int rv = opt.parse();
  if (rv)
    return 1;

  // decide where the output is going to appear
  setupOutputFile(opt);

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

  // close the output file
  if (opt.outputToFile())
    fclose(printout);

#ifdef WITH_MPI
  MPI_Finalize();
#endif

  return rv;
}

/**
 *  Decide where the output is going to appear.
 *
 *  Sets the global variable @c printout to decide where the output from
 *  OOPS is going to appear.
 *
 *  @param opt:
 *         Command line options.
 */
void setupOutputFile(OptionsOops &opt) {

  char filename[20];

  // redirect the output to a file
  if (opt.outputToFile()) {

#ifdef WITH_MPI
    sprintf(filename, "output%d.dat", MYID_PAR);
#else
    sprintf(filename, "output.dat");
#endif

    // open the file for writing
    printout = fopen(filename, "w");

#ifdef WITH_MPI
    fprintf(printout, "Output from processor %d.\n\n", MYID_PAR);
#endif
  }

  // print the output on the screen
  else
    printout = stdout;
}
