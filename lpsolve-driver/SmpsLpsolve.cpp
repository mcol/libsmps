/*
 *  SmpsLpsolve.cpp
 *
 *  SMPS interface to lpsolve.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include "SmpsLpsolve.h"
#include "lp_lib.h"


static int
createLP(lprec *lp, ProbData *prob);

static int
convertType(const int type);


/** Return status from solve() */
const char *solveStatus[10] = {
  "optimal",
  "suboptimal",
  "infeasible",
  "unbounded",
  "degenerate",
  "numerical failure",
  "user aborted",
  "timeout",
  "",
  "presolved"
};

/** Constructor */
SmpsLpsolve::SmpsLpsolve(string smpsFile) :
  smps(smpsFile) {}

/** Destructor */
SmpsLpsolve::~SmpsLpsolve() {

}

/** Read the smps files */
int SmpsLpsolve::read() {

  int rv = smps.read();
  if (rv)
    return rv;

  return rv;
}

/** Generate the deterministic equivalent */
ProbData* SmpsLpsolve::generateProblem() {

  return setupProblem(smps);
}

/** Solve the problem */
int SmpsLpsolve::solve(const OptionsLpsolve &opt) {

  int status = 0;

  lprec *lp = NULL;
  ProbData *prob = NULL;

  printf(" --------------- solve ---------------------\n");

  if (opt.writeMps() || opt.writeLp()) {
    smps.setBuildNames();
  }

  // generate the deterministic equivalent
  prob = generateProblem();
  if (!prob) {
    fprintf(stderr, "Failed to generate the deterministic equivalent.\n");
    status = 1;
    goto TERMINATE;
  }

  // create the problem
  lp = make_lp(prob->ttm + 1, prob->ttn); // also objective row counted here
  if (!lp) {
    fprintf(stderr, "Failed to create the problem.\n");
    status = 1;
    goto TERMINATE;
  }

  // fill the data of the LP problem
  status = createLP(lp, prob);
  if (status) {
    fprintf(stderr, "Failed to fill the data of the LP problem.\n");
    goto TERMINATE;
  }

  // write the deterministic equivalent in mps format
  if (opt.writeMps()) {
    status = write_freemps(lp, (char *) "smps.mps");
    if (status == 0) {
      fprintf(stderr, "Failed to write the mps file.\n");
      goto TERMINATE;
    }
  }

  // write the deterministic equivalent in lp format
  if (opt.writeLp()) {
    status = write_lp(lp, (char *) "smps.lp");
    if (status == 0) {
      fprintf(stderr, "Failed to write the lp file.\n");
      goto TERMINATE;
    }
  }

  // exit early if we don't have to solve the problem
  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // turn on presolve if it's explicitly required
  if (opt.usePresolve()) {
  set_presolve(lp, PRESOLVE_ROWS | PRESOLVE_COLS | PRESOLVE_LINDEP,
	       get_presolveloops(lp));
  }

  // solve the problem
  status = ::solve(lp);
  if (status) {
    fprintf(stderr, "Failed to optimize the LP problem (%s).\n",
	    solveStatus[status]);
    goto TERMINATE;
  }

  // get and print the solution
  status = getSolution(lp, opt);
  if (status) {
    fprintf(stderr, "Failed to retrieve the solution.\n");
    goto TERMINATE;
  }

 TERMINATE:

  // free up the problem as allocated by make_lp
  if (lp)
    delete_lp(lp);

  return status;
}

/** Fill the data of the LP problem */
int createLP(lprec *lp, ProbData *prob) {

  int status;

  int nRows = prob->ttm;
  int nCols = prob->ttn;
  int nNonz = prob->ttnz;

  double *obj = prob->obj;
  double *rhs = prob->rhs;
  double *blo = prob->blo;
  double *bup = prob->bup;
  char **rwnm = prob->rwnames;
  char **clnm = prob->clnames;

  // set the direction of optimization
  set_minim(lp);

  set_verbose(lp, 3);

  // shift the row numbers by 1, as row 0 should be left for the objective
  for (int i = 0; i < nNonz; ++i)
    ++prob->rwnmbs[i];

  // add the columns
  for (int col = 0; col < nCols; ++col) {
    int begCol = prob->clpnts[col];
    int endCol = prob->clpnts[col + 1];
    status = set_columnex(lp, col + 1, endCol - begCol,
			  &prob->acoeff[begCol], &prob->rwnmbs[begCol]);
    set_bounds(lp, col + 1, blo[col], bup[col]);
    if (!status)
      goto TERMINATE;
  }

  if (clnm) {
    for (int col = 0; col < nCols; ++col)
      set_col_name(lp, col + 1, clnm[col]);
  }

  // add the objective row
  {
    // lpsolve needs one extra element at the beginning of the objective
    double *obj1 = new double[nCols + 1];
    memcpy(&obj1[1], obj, nCols * sizeof(double));
    set_obj_fn(lp, obj1);
    delete[] obj1;
  }

  // set the signs of the constraints and the right-hand side
  for (int row = 1; row <= nRows; ++row) {
    set_constr_type(lp, row, convertType(prob->rws[row - 1]));
    set_rh(lp, row, rhs[row - 1]);
  }

  if (rwnm) {
    for (int row = 1; row <= nRows; ++row)
      set_row_name(lp, row, rwnm[row - 1]);
  }

 TERMINATE:

  // clean up
  freeProbData(prob);

  return !status;
}

/** Convert the type of the constraint to the lpsolve format */
int convertType(const int type)  {

  int newType;

  // convert the row codes to row types
  switch (type) {

    // equality
    case 1:
      newType = EQ;
      break;

    // >= inequality
    case 2:
      newType = GE;
      break;

    // <= inequality
    case 3:
      newType = LE;
      break;

    default:
      // we should not get here
      printf("Warning! Row code is %d.\n", type);
      newType = -1;
  }

  return newType;
}

/** Retrieve the solution information */
int SmpsLpsolve::getSolution(lprec *lp, const OptionsLpsolve &opt) {

  int status = 0, iters;
  double obj;

  int nRows = get_Norig_rows(lp);
  int nCols = get_Norig_columns(lp);
  double *primal = NULL, *dual = NULL, *slack = NULL, *rcosts = NULL;

  // retrieve the solution
  obj = get_objective(lp);

  // retrieve the number of iterations
  iters = get_total_iter(lp);

  // output the solution to the screen
  printf("\nObjective function value: %.8f (after %d iterations).\n",
	 obj, iters);

  if (opt.printSolution()) {

    // this contains the value of the objective function [0], of the
    // constraints [1, nRows], and of the variables [nRows + 1, nRows + nCols]
    double *psol = new double[1 + nRows + nCols];
    get_primal_solution(lp, psol);

    // this contains a null value [0], the duals of the constraints [1, nRows],
    // and the duals of the variables [nRows + 1, nRows + nCols]
    double *dsol = new double[1 + nRows + nCols];
    get_dual_solution(lp, dsol);

    // compute the slack variables
    for (int row = 1; row <= nRows; ++row)
      psol[row] = get_rh(lp, row) - psol[row];

    primal = &psol[nRows + 1];
    dual   = &dsol[1];
    slack  = &psol[1];
    rcosts = &dsol[nRows + 1];

    printSolution(smps.getRootNode(), primal, dual, slack, rcosts);
    delete[] psol;
    delete[] dsol;
  }

  return status;
}

/** Constructor */
OptionsLpsolve::OptionsLpsolve(const int argc, const char *argv[]) :
  Options(argc, argv),
  _writeLp(0),
  _usePresolve(0) {
}

/** Parse the command line options */
int OptionsLpsolve::parse() {

  // add the specialised options
  Options::addOption("-l", "write the deterministic equivalent in LP format",
		     &_writeLp);
  Options::addOption("-p", "turn on the presolve", &_usePresolve);

  // parse the common options
  int rv = Options::parse();

  return rv;
}
