/*
 *  SmpsGlpk.cpp
 *
 *  SMPS interface to Glpk.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include "SmpsGlpk.h"


static int
createLP(glp_prob *lp, ProbData *prob);

static int
convertColPointers(int *colPointers, int *colNumbers, const int nCols);


/** Constructor */
SmpsGlpk::SmpsGlpk(string smpsFile) :
  smps(smpsFile) {}

/** Destructor */
SmpsGlpk::~SmpsGlpk() {

}

/** Read the smps files */
int SmpsGlpk::read() {

  int rv = smps.read();
  if (rv)
    return rv;

  return rv;
}

/** Generate the deterministic equivalent */
ProbData* SmpsGlpk::generateProblem() {

  return setupProblem(smps);
}

/** Solve the problem */
int SmpsGlpk::solve(const OptionsGlpk &opt) {

  int status = 0;

  glp_prob *lp = NULL;
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
  lp = glp_create_prob();
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
    status = glp_write_mps(lp, GLP_MPS_FILE, NULL, "smps.mps");
    if (status) {
      fprintf(stderr, "Failed to write the mps file.\n");
      goto TERMINATE;
    }
  }

  // write the deterministic equivalent in lp format
  if (opt.writeLp()) {
    status = glp_write_lp(lp, NULL, "smps.lp");
    if (status) {
      fprintf(stderr, "Failed to write the lp file.\n");
      goto TERMINATE;
    }
  }

  // exit early if we don't have to solve the problem
  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // solve the problem
  if (opt.useBarrier())
    status = lpx_interior(lp);
  else
    status = glp_simplex(lp, NULL);

  if (status && status != LPX_E_OK) {
    fprintf(stderr, "Failed to optimize the LP problem.\n");
    goto TERMINATE;
  }

  // get and print the solution
  status = getSolution(lp, opt);
  if (status) {
    fprintf(stderr, "Failed to retrieve the solution.\n");
    goto TERMINATE;
  }

 TERMINATE:

  // free up the problem as allocated by glp_create_prob
  if (lp)
    glp_delete_prob(lp);

  return status;
}

/** Fill the data of the LP problem */
int createLP(glp_prob *lp, ProbData *prob) {

  int nRows = prob->ttm;
  int nCols = prob->ttn;
  int nNonz = prob->ttnz;

  double *obj = prob->obj;
  double *rhs = prob->rhs;
  double *blo = prob->blo;
  double *bup = prob->bup;
  char **rwnm = prob->rwnames;
  char **clnm = prob->clnames;

  // array of column numbers
  int *colnmbs = new int[nNonz];

  // set the direction of optimization
  glp_set_obj_dir(lp, GLP_MIN);

  // add rows and columns
  glp_add_rows(lp, nRows);
  glp_add_cols(lp, nCols);

  // convert clpnts to colnmbs
  convertColPointers(prob->clpnts, colnmbs, nCols);

  // row and column numbers should start from 1
  for (int i = 0; i < nNonz; ++i) {
    ++prob->rwnmbs[i];
    ++colnmbs[i];
  }

  // load the matrix data
  // since glpk expects the arrays to have dimension nNonz + 1 and the
  // first coefficient is ignored, we pass the address of the arrays - 1
  glp_load_matrix(lp, nNonz, prob->rwnmbs - 1, colnmbs - 1, prob->acoeff - 1);

  // set up the sense of inequalities and right-hand sides
  for (int row = 0; row < nRows; ++row) {

    const int type = prob->rws[row];

    switch (type) {

      // equality
      case 1:
	glp_set_row_bnds(lp, row + 1, GLP_FX, rhs[row], rhs[row]);
	break;

      // >= inequality
      case 2:
	glp_set_row_bnds(lp, row + 1, GLP_LO, rhs[row], 0.0);
	break;

      // <= inequality
      case 3:
	glp_set_row_bnds(lp, row + 1, GLP_UP, 0.0, rhs[row]);
	break;

      default:
	// we should not get here
	printf("Warning! Row code is %d.\n", type);
    }
  }

  // set up the objective and the bounds
  for (int col = 0; col < nCols; ++col) {

    glp_set_obj_coef(lp, col + 1, obj[col]);

    if (blo[col] < -1.e20 && bup[col] > 1.e20) {
      printf("%d is free var\n", col);
      glp_set_col_bnds(lp, col + 1, GLP_FR, blo[col], bup[col]);
    }
    else if (bup[col] < 1.e20) {
      printf("%d is upbound var\n", col);
      glp_set_col_bnds(lp, col + 1, GLP_DB, 0.0, bup[col]);
    }
    else
      glp_set_col_bnds(lp, col + 1, GLP_LO, blo[col], 0.0);
  }

  // assign the row names
  if (rwnm) {
    for (int row = 0; row < nRows; ++row)
      glp_set_row_name(lp, row + 1, rwnm[row]);
  }

  // assign the column names
  if (clnm) {
    for (int col = 0; col < nCols; ++col)
      glp_set_col_name(lp, col + 1, clnm[col]);
  }

  // clean up
  delete[] colnmbs;
  freeProbData(prob);

  return 0;
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

/** Retrieve the solution information */
int SmpsGlpk::getSolution(glp_prob *lp, const OptionsGlpk &opt) {

  double obj;

  // retrieve the solution
  if (opt.useBarrier())
    obj = glp_ipt_obj_val(lp);
  else
    obj = glp_get_obj_val(lp);

  // output the solution to the screen
  printf("\nObjective function value: %.8f.\n", obj);

  if (opt.printSolution()) {

    const int nRows = glp_get_num_rows(lp);
    const int nCols = glp_get_num_cols(lp);
    double *primal = new double[nCols];
    double *dual   = new double[nRows];
    double *slack  = new double[nRows];
    double *rcosts = new double[nCols];
    double rhs;

    for (int row = 0; row < nRows; ++row) {
      dual[row]  = glp_get_row_dual(lp, row + 1);

      // compute the slack
      int rwtype = glp_get_row_type(lp, row + 1);
      if (rwtype == GLP_LO)
	rhs = glp_get_row_lb(lp, row + 1);
      else if (rwtype == GLP_UP)
	rhs = glp_get_row_ub(lp, row + 1);
      printf("rhs]%d]: %f\n", row, rhs);
      slack[row] = rhs - glp_get_row_prim(lp, row + 1);
    }

    for (int col = 0; col < nCols; ++col) {
      primal[col] = glp_get_col_prim(lp, col + 1);
      rcosts[col] = glp_get_col_dual(lp, col + 1);
    }


    // this contains the value of the objective function [0], of the
    // constraints [1, nRows], and of the variables [nRows + 1, nRows + nCols]
    //    double *psol = new double[1 + nRows + nCols];
    //    get_primal_solution(lp, psol);

    // this contains a null value [0], the duals of the constraints [1, nRows],
    // and the duals of the variables [nRows + 1, nRows + nCols]
    //    double *dsol = new double[1 + nRows + nCols];
    //    get_dual_solution(lp, dsol);

    // compute the slack variables
    //    for (int row = 1; row <= nRows; ++row)
    //      psol[row] = get_rh(lp, row) - psol[row];

    printSolution(smps.getRootNode(), primal, dual, slack, rcosts);
    delete[] primal;
    delete[] dual;
    delete[] slack;
    delete[] rcosts;
  }

  return 0;
}

/** Constructor */
OptionsGlpk::OptionsGlpk(const int argc, const char *argv[]) :
  Options(argc, argv),
  _writeLp(0),
  _useBarrier(0) {
}

/** Parse the command line options */
int OptionsGlpk::parse() {

  // add the specialised options
  Options::addOption("-l", "write the deterministic equivalent in LP format",
		     &_writeLp);
  Options::addOption("-b", "use a barrier solver", &_useBarrier);

  // parse the common options
  int rv = Options::parse();

  return rv;
}
