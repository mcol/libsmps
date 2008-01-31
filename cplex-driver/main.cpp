/*
 *  main.cpp
 *
 *  Driver for the SMPS interface to Cplex.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdlib.h>
#include "SmpsCplex.h"


static int
createLP(CPXENVptr env, CPXLPptr lp, ProbData *prob);

static int
setRowType(int *rowCode, char *rowType, const int nRows);

static int
convertColPointers(int *colPointers, int *colNumbers, const int nCols);


/** Driver routine for the SMPS interface to Cplex */
int main(const int argc, const char *argv[]) {

  // parse the command line options
  OptionsCplex opt(argc, argv);
  int rv = opt.parse();
  if (rv)
    return 1;

  // create an object for the problem
  SmpsCplex data(opt.getSmpsFile());

  // read the smps files
  rv = data.read();
  if (rv)
    goto TERMINATE;

  // solve the problem
  rv = data.solve(opt);

 TERMINATE:

  return rv;
}

/** Solve the problem */
int SmpsCplex::solve(const OptionsCplex &opt) {

  int status = 0;

  CPXENVptr env;
  CPXLPptr lp;

  printf(" --------------- cplexSolve ----------------\n");

  // initialize the CPLEX environment
  env = CPXopenCPLEX(&status);
  if (!env) {

    char errmsg[1024];

    CPXgeterrorstring(env, status, errmsg);
    fprintf(stderr, "Could not open the CPLEX environment.\n%s", errmsg);

    goto TERMINATE;
  }

  // turn on output to the screen
  status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
  if (status) {
    fprintf(stderr, "Failed to turn on screen indicator, error %d.\n", status);
    goto TERMINATE;
  }

  // turn off presolve unless it's explicitly required
  if (!opt.usePresolve()) {
    status = CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);
    if (status) {
      fprintf(stderr, "Failed to turn off presolve, error %d.\n", status);
      goto TERMINATE;
    }
  }

  // create the problem
  lp = CPXcreateprob(env, &status, "lp");
  if (!lp) {
    fprintf(stderr, "Failed to create the problem.\n");
    goto TERMINATE;
  }

  // fill the data of the LP problem
  ProbData *prob = generateProblem();
  status = createLP(env, lp, prob);
  if (status) {
    fprintf(stderr, "Failed to fill the data of the LP problem.\n");
    goto TERMINATE;
  }

  // write the deterministic equivalent in mps format
  if (opt.writeMps()) {
    status = CPXwriteprob(env, lp, "smps.mps", NULL);
    if (status) {
      fprintf(stderr, "Failed to write the mps file.\n");
      return 1;
    }
  }

  // leave early if we don't have to solve the problem
  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // solve the problem
  if (opt.useBarrier())
    status = CPXbaropt(env, lp);
  else
    status = CPXlpopt(env, lp);

  if (status) {
    fprintf(stderr, "Failed to optimize the LP problem.\n");
    goto TERMINATE;
  }

  // get and print the solution
  status = getSolution(env, lp, opt);
  if (status) {
    fprintf(stderr, "Failed to obtain the solution.\n");
    goto TERMINATE;
  }

 TERMINATE:

  // free up the problem as allocated by CPXcreateprob, if necessary
  if (lp) {
    status = CPXfreeprob(env, &lp);
    if (status)
      fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
  }

  // free up the CPLEX environment, if necessary
  if (env) {
    status = CPXcloseCPLEX(&env);

    if (status) {

      char errmsg[1024];

      CPXgeterrorstring(env, status, errmsg);
      fprintf(stderr, "Could not close the CPLEX environment.\n%s", errmsg);
    }
  }

  return status;
}

/** Fill the data of the LP problem */
int createLP(CPXENVptr env, CPXLPptr lp, ProbData *prob) {

  int status;

  int nRows = prob->ttm;
  int nCols = prob->ttn;
  int nNonz = prob->ttnz;

  double *obj = prob->obj;
  double *rhs = prob->rhs;
  double *blo = prob->blo;
  double *bup = prob->bup;

  // array of column numbers
  int *colnmbs = new int[nNonz];

  // array for the direction of the inequalities
  char *rowType = new char[nRows];

  // set the direction of optimization
  CPXchgobjsen(env, lp, CPX_MIN);

  // convert clpnts to colnmbs
  convertColPointers(prob->clpnts, colnmbs, nCols);

  // decide the signs of the constraints
  status = setRowType(prob->rws, rowType, nRows);
  if (status)
    goto TERMINATE;

  // create the new rows
  status = CPXnewrows(env, lp, nRows, rhs, rowType, NULL, NULL);
  if (status)
    goto TERMINATE;

  // add the columns
  status = CPXnewcols(env, lp, nCols, obj, blo, bup, NULL, NULL);
  if (status)
    goto TERMINATE;

  // create the list of coefficients
  status = CPXchgcoeflist(env, lp, nNonz, prob->rwnmbs, colnmbs, prob->acoeff);
  if (status)
    goto TERMINATE;

 TERMINATE:

  // clean up
  delete[] colnmbs;
  delete[] rowType;
  freeProbData(prob);

  return status;
}

/** Convert the column pointers */
int convertColPointers(int *colPointers, int *colNumbers, const int nCols) {

  // for all columns
  for (int col = 0; col < nCols; ++col) {

    for (int i = colPointers[col]; i < colPointers[col + 1]; ++i)
      colNumbers[i] = col;
  }

  return 0;
}

/** Set the type of the constraints */
int setRowType(int *rowCode, char *rowType, const int nRows)  {

  for (int i = 0; i < nRows; ++i) {

    // convert the row codes to row types
    switch (rowCode[i]) {

      // equality
      case 1:
	rowType[i] = 'E';
	break;

      // >= inequality
      case 2:
	rowType[i] = 'G';
	break;

      // <= inequality
      case 3:
	rowType[i] = 'L';
	break;

      default:
	// we should not get here
	printf("Warning! Row code is %d.\n", rowCode[i]);
	return 1;
    }
  }

  return 0;
}

/** Retrieve the solution information */
int SmpsCplex::getSolution(CPXENVptr env, CPXLPptr lp,
			   const OptionsCplex &opt) {

  int status, solstat, iters;
  double obj;

  int numRows = CPXgetnumrows(env, lp);
  int numCols = CPXgetnumcols(env, lp);

  double *primal = new double[numCols];
  double *dual   = new double[numRows];
  double *slack  = new double[numRows];
  double *rcosts = new double[numCols];

  if (primal == NULL || dual   == NULL ||
      slack  == NULL || rcosts == NULL) {
    status = CPXERR_NO_MEMORY;
    fprintf(stderr, "Could not allocate memory for the solution.\n");
    goto TERMINATE;
  }

  // retrieve the solution
  status = CPXsolution(env, lp, &solstat, &obj, primal, dual, slack, rcosts);
  if (status)
    goto TERMINATE;

  // retrieve the number of iterations
  if (CPXgetmethod(env, lp) == CPX_ALG_BARRIER)
    iters = CPXgetbaritcnt(env, lp);
  else
    iters = CPXgetitcnt(env, lp);

  // output the solution to the screen
  printf("\nObjective function value: %.8f (after %d iterations).\n",
	 obj, iters);

  if (opt.printSolution()) {
    const NodeInfo *info = getNodeInfo();
    printSolution(info, primal, dual, slack, rcosts);
    delete info;
  }

 TERMINATE:

  // clean up
  delete[] primal;
  delete[] dual;
  delete[] slack;
  delete[] rcosts;

  return status;
}
