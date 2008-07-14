/*
 *  main.cpp
 *
 *  Driver for the SMPS interface to Hopdm.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "SmpsHopdm.h"


/** Driver routine for the SMPS interface to Hopdm */
int main (const int argc, const char *argv[]) {

  // parse the comman line options
  OptionsHopdm opt(argc, argv);
  int rv = opt.parse();
  if (rv)
    return 1;

  // create an object for the problem
  SmpsHopdm data(opt.smpsFile());

  // read the smps files
  rv = data.read();
  if (rv)
    goto TERMINATE;

#ifdef WITH_WARMSTART
  // warmstart case
  if (opt.useWarmstart()) {

    // create a reduced tree
    if (opt.useReduction())
      rv = data.reduceTree(opt.useReduction());

    if (rv)
      goto LEAVE_WARMSTART;

    // solve the reduced tree
    rv = data.solveReduced(opt);
    if (rv)
      goto TERMINATE;
  }

 LEAVE_WARMSTART:
#endif

  // solve the problem
  rv = data.solve(opt);

 TERMINATE:

  return rv;
}

#if 0

  NodeInfo *info = NULL;

  Ws **warmStart = (Ws **) calloc(1, sizeof(Ws *));

  // set the mode of operation
  int mode = options->warmstart ? 1 : 0;

  // generate and solve a reduced problem to obtain a warmstarting point
  if (options->warmstart) {

    int failed;

    SmpsData* rData = (SmpsData *) calloc(1, sizeof(SmpsData));

    // choose the scenarios to appear in the reduced tree
    int *chosen = chooseScenarios(data->tree, 6);

    // generate the reduced tree
    SmpsTree *rTree = generateReducedTree(data, chosen);

    rData->core = data->core;
    rData->time = data->time;
    rData->scen = data->scen;
    rData->tree = rTree;

    // generate the reduced problem
    ProbData *rProb = setupProblem(rData);

    // solve the reduced problem
    failed = hopdmSolve(rProb, warmStart, mode, NULL);
    mode = 2;

    // build up the warm-start iterate if an optimal solution has been found
    if (!failed)
      setupWarmStart(data->tree, rTree, warmStart, chosen);

    // free the reduced tree
    freeSmpsTree(rTree);

    // clean up
    free(chosen);
    free(rData);

    // clean up for the early termination
    if (failed) {
      freeSmpsData(data);
      goto TERMINATE;
    }
  }

  // store the information necessary to print the solution node by node
  if (options->printSoln)
    info = newNodeInfo(data->tree);

  // generate the problem
  ProbData *prob = setupProblem(data);

  // write the deterministic equivalent in mps format
  if (options->writeMps)
    writeMps("smps.mps", prob);

  // free the smps data structure
  freeSmpsData(data);

  // solve the problem
  hopdmSolve(prob, warmStart, mode, info);

 TERMINATE:

  // clean up
  freeOptions(options);
  free(warmStart);
  if (info)
    freeNodeInfo(info);

  return 0;
}

#endif // 0
