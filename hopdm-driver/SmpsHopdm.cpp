/*
 *  SmpsHopdm.cpp
 *
 *  SMPS interface to Hopdm.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "SmpsHopdm.h"

static int
hopdm(Lp *lp, Ws *warmStart);

static int
addObjectiveRow(ProbData *prob);

static Lp*
defineProblem(ProbData *prob, const int mode);

#if 0
static int wspdm(Lp *lp, Ws *warmStart);
static char* writeRowNames(const int nRows);
#endif


/** Constructor */
SmpsHopdm::SmpsHopdm(string smpsFile) :
  smps(smpsFile),
  rTree() {
}

/** Read the smps files */
int SmpsHopdm::read() {

  int rv = smps.read();
  if (rv)
    return rv;

  return rv;
}

/** Generate the deterministic equivalent */
ProbData* SmpsHopdm::generateProblem(const Node *root) {

  return setupProblem(smps, root);
}

/** Solve the problem */
int SmpsHopdm::solve(const OptionsHopdm &opt) {

  int rv = 0;
  Lp *lp = NULL;
  ProbData *prob = NULL;

  printf(" --------------- solve ---------------------\n");

  if (opt.writeMps()) {
    // disabled because hopdm wants the objective row to be included in the
    // matrix, which causes an off by 1 error that appears in valgrind during
    // freeProbData() (getTotRows() and prob->ttm differ by 1 element)
    printf("The writing of the MPS file is disabled.\n");
  }

  // generate the deterministic equivalent
  prob = generateProblem(smps.getRootNode());
  if (!prob) {
    fprintf(stderr, "Failed to generate the deterministic equivalent.\n");
    rv = 1;
    goto TERMINATE;
  }

  // add the objective row to the matrix
  rv = addObjectiveRow(prob);
  if (rv)
    goto TERMINATE;

  // define the problem
  lp = defineProblem(prob, 0);
  if (!lp) {
    rv = 1;
    goto TERMINATE;
  }

  // exit early if we don't have to solve the problem
  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // solve the problem
  rv = hopdm(lp, NULL);
  if (rv) {
    fprintf(stderr, "Failed to optimize the LP problem.\n");
    goto TERMINATE;
  }

  // get and print the solution
  rv = getSolution(lp, opt);
  if (rv) {
    fprintf(stderr, "Failed to retrieve the solution.\n");
    goto TERMINATE;
  }

 TERMINATE:

  // clean up
  free_Lp(lp);

  return rv;
}

#if 0
/** Solve a reduced problem */
int SmpsHopdm::solveReduced(const OptionsHopdm &opt) {

  int rv = 0;
  Lp *lp = NULL;
  ProbData *prob = NULL;

  printf(" --------------- solveReduced --------------\n");

  // generate a reduced problem
  prob = generateProblem(rTree.getRootNode());
  if (!prob) {
    printf("Failed to generate the deterministic equivalent.\n");
    rv = 1;
    goto TERMINATE;
  }

  // add the objective row to the matrix
  rv = addObjectiveRow(prob);
  if (rv)
    goto TERMINATE;

  // define the problem
  lp = defineProblem(prob, 0);
  if (!lp) {
    rv = 1;
    goto TERMINATE;
  }

  // exit early if we don't have to solve the problem
  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // extract and store the solution
  storeSolution(pdProb, prob);

  // generate a warmstart point for the complete problem
  setupWarmStart(prob);

 TERMINATE:

  // clean up
  free_Lp(lp);

  return rv;
}

/* -------------------------------------------------------------------------
 * hopdmSolve
 *
 * The mode determines what algorithm to use:
 *   0: solve the problem from scratch
 *   1: find a warmstart iterate
 *   2: solve the problem from the warmstart iterate
 * ------------------------------------------------------------------------- */
int hopdmSolve(ProbData *prob, Ws **iterate, const int mode) {

  int code;
  Lp *lp;

  printf(" --------------- hopdmSolve ----------------\n");

  // add the objective row to the matrix
  code = addObjectiveRow(prob);

  // allocate space for the warmstart iterate
  if (!*iterate)
    *iterate = (Ws *) calloc(1, sizeof(Ws));

  // define the problem
  lp = defineProblem(prob, mode);
  if (!lp)
    return 1;

  // solve the problem
  if (mode == 1)
    // find the warmstart iterate
    code = wspdm(lp, *iterate);
  else
    // find the optimal solution
    code = hopdm(lp, *iterate);

  // clean up
  free_Lp(lp);

  return code;
}
#endif

/* -------------------------------------------------------------------------
 * hopdm
 * ------------------------------------------------------------------------- */
int hopdm(Lp *lp, Ws *warmStart) {

  int status;
  int restart;
  long iexit = 0;        // get optimal solution (no matter of centering)
  long iters;

  double opttol = 5.0e-8;
  double errb, erru, errc, asmall, alarge;

  if (!warmStart) {
    warmStart = new Ws;
    memset(warmStart, 0, sizeof(Ws));
  }

  // solve the problem from scratch
  if (!warmStart->xsave) {
    restart = 0;
    status = solve_Lp(lp, restart, warmStart, &opttol, &iexit,
		    &errb, &erru, &errc, &asmall, &alarge, &iters);
  }

  // solve the problem from a warmstart iterate
  else {
    restart = 1;
    printf("Solving the complete problem.\n");
    status = wspdm_Lp(lp, restart, warmStart, &opttol, &iexit,
		      &errb, &erru, &errc, &asmall, &alarge, &iters, NULL);
  }

  if (status) {
    goto TERMINATE;
  }

  // store the number of iterations in the Lp structure
  lp->itslv = iters;

 TERMINATE:

  // clean up
  if (warmStart->xsave)
    free_Lp_save(warmStart);
  delete warmStart;

  return status;
}

#if 0
/** */
int wspdm(Lp *lp, Ws *warmStart) {

  int status;
  int restart = 0;        // cold start
  long iexit  = 1;        // get a centered solution
  long iters;

  double opttol = 5.0e-1;
  double errb, erru, errc, asmall, alarge;

  // solve the problem to find a warmstart iterate
  printf("Solving the reduced problem.\n");
  status = wspdm_Lp(lp, restart, warmStart, &opttol, &iexit,
		    &errb, &erru, &errc, &asmall, &alarge, &iters, NULL);

  if (status) {
    fprintf(stderr, "Failed to optimize the LP problem, error %d.\n", status);
    return status;
  }

  // optimal solution found
  //  getSolution(lp, opt);  // moved inside solve()

  // store the optimal solution
  save_Lp(lp, warmStart);

  return 0;
}
#endif

/** Hopdm expects the objective row as the last row of the matrix */
int addObjectiveRow(ProbData *prob) {

  int col, endCol;
  int oldIndex = 0, newIndex = 0;

  int nRows  = prob->ttm;
  int nCols  = prob->ttn;
  int nNonz  = prob->ttnz;

  // put the objective in the last row of the matrix
  int objRow = prob->ttm;

  double *obj    = prob->obj;
  double *acoeff = prob->acoeff;
  int    *rwnmbs = prob->rwnmbs;
  int    *clpnts = prob->clpnts;

  // allocate enough additional space for a full objective
  double *newacoeff = new double[nNonz + nCols];
  int    *newrwnmbs = new int[nNonz + nCols];

  // reallocate the space for rws
  int *newrws = new int[nRows + 1];
  memcpy(newrws, prob->rws, nRows *  sizeof(int));
  delete[] prob->rws;
  prob->rws = newrws;

  // reallocate the space for rhs
  double *newrhs = new double[nRows + 1];
  memcpy(newrhs, prob->rhs, nRows *  sizeof(double));
  delete[] prob->rhs;
  prob->rhs = newrhs;

  // for all columns
  for (col = 0; col < nCols; col++) {

    assert(oldIndex < nNonz);

    // set the new column pointer
    clpnts[col] = newIndex;

    // index of the last element in the current column
    endCol = clpnts[col + 1];

    // copy all elements in the current column
    while (oldIndex < endCol) {

      newacoeff[newIndex] = acoeff[oldIndex];
      newrwnmbs[newIndex] = rwnmbs[oldIndex];
      oldIndex++;
      newIndex++;
    }

    // copy the objective, if nonzero
    if (obj[col]) {

      newacoeff[newIndex] = obj[col];
      newrwnmbs[newIndex] = objRow;
      newIndex++;
    }
  }

  // set the last column pointer
  clpnts[col] = newIndex;

  assert(newIndex <= nNonz + nCols);
  assert(clpnts[nCols] == newIndex);

  // clean up the old arrays
  delete[] acoeff;
  delete[] rwnmbs;

  // point to the new arrays
  prob->acoeff = newacoeff;
  prob->rwnmbs = newrwnmbs;

  // set the information about the objective row
  prob->rws[objRow] = 4;
  prob->rhs[objRow] = 0.0;
  prob->irobj = objRow;

  // update the sizes of the problem
  prob->ttm += 1;
  prob->ttnz = newIndex;

  return 0;
}

/** */
Lp* defineProblem(ProbData *prob, const int mode) {

  int i, status;
  int iolog    = 77;
  int presolve = 0;

  double mult = 1.0;
  double big  = 1.e30;
  double objcon = 0.0;  // objective constant

  int m   = prob->ttm;
  int n   = prob->ttn;
  int nza = prob->ttnz;
  int obj = prob->irobj;

  int max_rows = m + 1;
  int max_cols = n + 1;
  int max_nza  = nza + 1;
  int max_nzq  = 1;

  long *clpnts = (long *) prob->clpnts;
  long *rwnmbs = (long *) prob->rwnmbs;
  long *rwstat = (long *) prob->rws;
  double *acf = prob->acoeff;
  double *rhs = prob->rhs;
  double *blo = prob->blo;
  double *bup = prob->bup;
  double *ranges = (double *) calloc(m, sizeof(double));

  MPS *P;
  Lp  *lpProblem = NULL;

  // set the ranges
  for (i = 0; i < m; i++)
    ranges[i] = big;

  // fortran numbering for obj, clpnts and rwnmbs
  obj++;
  for (i = 0; i < n + 1; i++)
    clpnts[i]++;
  for (i = 0; i < nza; i++)
    rwnmbs[i]++;

  // switch the presolve off if we are dealing with warmstart
  if (mode)
    presolve = 0;

  P = init_mps(max_rows, max_cols, max_nza, max_nzq,
	       presolve, iolog, mult, big);

  status = def_MPS_dimensions(P, m, n);
  if (status) {
    fprintf(stderr, "Failed to define the MPS dimensions.\n");
    goto TERMINATE;
  }

  status = define_A(P, nza, obj, clpnts, rwnmbs, acf);
  if (status) {
    fprintf(stderr, "Failed to define matrix A.\n");
    goto TERMINATE;
  }

  status = def_MPS_data(P, objcon, rhs, ranges, rwstat, blo, bup);
  if (status) {
    fprintf(stderr, "Failed to define the MPS data.\n");
    goto TERMINATE;
  }

  status = define_names(P, NULL, NULL);
  if (status) {
    fprintf(stderr, "Failed to define the MPS names.\n");
    goto TERMINATE;
  }

  lpProblem = mps2hopdm(P);
  if (!lpProblem) {
    fprintf(stderr, "Failed to create the LP problem.\n");
    goto TERMINATE;
  }

  // preprocessing
  preproc_Lp(lpProblem);

 TERMINATE:

  // clean up
  free(ranges);
  free_mps(P);
  freeProbData(prob);

  return lpProblem;
}

/** Retrieve the solution information */
int SmpsHopdm::getSolution(Lp *lp, const OptionsHopdm &opt) {

  print_Lp(lp, (char *) "solution.mps", (char *) "LP", 0, NULL);

  printf("\nObjective function value: %.8f (after %ld iterations).\n",
	 lp->objective, lp->itslv);

  if (!opt.printSolution())
    return 0;

  // allocate space for the solution vectors
  double *primal = (double *) calloc(lp->ncopy, sizeof(double));
  double *dual   = (double *) calloc(lp->mcopy, sizeof(double));
  double *slack  = (double *) calloc(lp->mcopy, sizeof(double));
  double *rcosts = (double *) calloc(lp->ncopy, sizeof(double));
  long   status;
  double objective;

  get_opt_solution(lp, &status, &objective, primal, dual, slack, rcosts);
  printSolution(smps.getRootNode(), primal, dual, slack, rcosts);

  // clean up
  free(primal);
  free(dual);
  free(slack);
  free(rcosts);

  return 0;
}

#if 0
/**Create a warmstart iterate from the reduced tree solution. */
void setupWarmStart(const SmpsTree *cTree, const SmpsTree *rTree,
		    Ws **iterate, const int *chosen) {

  int node;
  int nRows = cTree->f_rw_nd[cTree->nNodes] + 1;   // objective included
  int nCols = cTree->f_cl_nd[cTree->nNodes];
  Ws *rIterate = *iterate;
  Ws *new = (Ws *) calloc(1, sizeof(Ws));

  printf(" --------------- setupWarmStart ------------\n");

  // set the dimensions of the complete iterate
  new->m  = nRows - 1;
  new->m0 = nRows;
  new->n  = 1;
  new->n_struct = nCols;
  new->bsave = rIterate->bsave;

  // allocate space for the vectors in the complete iterate
  new->xsave  = (double *) calloc(nCols, sizeof(double));
  new->zsave  = (double *) calloc(nCols, sizeof(double));
  new->ssave  = (double *) calloc(nCols, sizeof(double));
  new->wsave  = (double *) calloc(nCols, sizeof(double));
  new->ysave  = (double *) calloc(nRows, sizeof(double));
  new->rwsave = writeRowNames(nRows);

  for (node = 0; node <= cTree->nNodes; node++)
    printf("%2d (%2d)  %3d\n", node, cTree->order[node], cTree->f_rw_nd[node]);

  // fill up the complete arrays
  for (node = 0; node < cTree->nNodes; node++) {

    int i, cIndex, rIndex, nElems;
    int cNode = cTree->order[node];
    int rNode = chosen[findNode(cTree, chosen, cNode)] - 1;

    // copy the x, z, s and w iterates
    cIndex = cTree->f_cl_nd[cNode];
    rIndex = rTree->f_cl_nd[rNode];

    nElems = cTree->f_cl_nd[cNode + 1] - cIndex;
    /*
    if (node + 1 == cTree->nNodes)
      nElems = cTree->f_cl_nd[node + 1] - cIndex;
    else
      nElems = cTree->f_cl_nd[cTree->order[node + 1]] - cIndex;
    */
    assert(nElems > 0);

    /*
    printf("cNode: %2d, cInded: %2d, rNode: %2d, rIndex: %2d, nElems: %2d\n",
	   cNode + 1, cIndex, rNode + 1, rIndex, nElems);
    */

    for (i = 0; i < nElems; i++) {

      assert(cIndex + i < nCols);
      assert(rIndex + i < rTree->f_cl_nd[rTree->nNodes]);
      new->xsave[cIndex + i] = rIterate->xsave[rIndex + i];
      new->zsave[cIndex + i] = rIterate->zsave[rIndex + i];
      new->ssave[cIndex + i] = rIterate->ssave[rIndex + i];
      new->wsave[cIndex + i] = rIterate->wsave[rIndex + i];
    }

    // copy the y iterate
    cIndex = cTree->f_rw_nd[cNode];
    rIndex = rTree->f_rw_nd[rNode];
    nElems = cTree->f_rw_nd[cNode + 1] - cIndex;

    /*
    if (node + 1 == cTree->nNodes)
      nElems = cTree->f_rw_nd[node + 1] - cIndex;
    else
      nElems = cTree->f_rw_nd[cTree->order[node + 1]] - cIndex;
    */
    assert(nElems > 0);

#ifdef DEBUG
    printf("%d: cNode: %d (rNode: %d) nElems: %d - cIndex: %d, rIndex: %d\n",
	   node, cNode, rNode, nElems, cIndex, rIndex);
#endif

    for (i = 0; i < nElems; i++) {
      new->ysave[cIndex + i] = rIterate->ysave[rIndex + i] * cTree->probnd[cNode];
    }
  }

#ifdef DEBUG
  int i;
  //  for (i = 0; i < nCols; i++)
  //    printf("%2d: %12.7f  %12.7f\n", i, new->xsave[i], new->zsave[i]);

  for (i = 0; i < nRows - 1; i++)
    printf("%2d: %12.7f\n", i, new->ysave[i]);
#endif

  // clean up the reduced iterate
  free_Lp_save(rIterate);
  free(rIterate);

  // point to the complete iterate
  *iterate = new;
}

/** */
char* writeRowNames(const int nRows) {

  int row;
  char *name = (char *) calloc(nRows*8 + 1, 1);
  char *temp = new char[8 + 1];
  char *pos  = name;

  for (row = 1; row <= nRows; row++) {
    sprintf(temp, "R%-7d", row);
    memcpy(pos, temp, 8*sizeof(char));
    pos += 8;
  }

  // clean up
  delete[] temp;

  return name;
}
#endif

/** Constructor */
OptionsHopdm::OptionsHopdm(const int argc, const char *argv[]) :
  Options(argc, argv),
  _useReduction(0),
  _useAggregation(0) {
}

/** Parse the command line options */
int OptionsHopdm::parse() {

  // add the specialised options
  Options::addOption("-w", "use a warmstart strategy with scenario reduction",
                     &_useReduction, true);
  Options::addOption("-a", "use a warmstart strategy with stage aggregation",
                     &_useAggregation, true);

  // parse the common options
  int rv = Options::parse();

  return rv;
}
