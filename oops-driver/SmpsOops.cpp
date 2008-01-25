/*
 *  SmpsOops.cpp
 *
 *  SMPS interface to Oops.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include "SmpsOops.h"
#include "oops/oopstime.h"
#include "oops/WriteMps.h"


double tt_start, tt_end;


/** Constructor */
SmpsOops::SmpsOops(string smpsFile, const int lev) :
  smps(smpsFile),
  level(lev),
  nBlocks(0),
  block(NULL),
  order(NULL),
  revorder(NULL) {
}

/** Destructor */
SmpsOops::~SmpsOops() {

  delete[] block;
  delete[] order;
  delete[] revorder;
}

/** Read the smps files */
int SmpsOops::read() {

  int rv = smps.read();
  if (rv)
    return rv;

  // reorder the nodes according to the level
  rv = orderNodes();

  return rv;
}

/**
 *  Call OOPS to solve the problem.
 *
 *  @param opt:
 *         Command line options.
 *  @param hopdm_options:
 *         Pointer to the options for the solver.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::solve(const OptionsOops &opt, HopdmOptions *hopdm_options) {

  // generate the problem
  SmpsReturn *prob = generateSmps();

  // check the validity of the arguments
  if (!prob) {
    cout << "No problem defined." << endl;
    return 1;
  }

  // setup the primal-dual problem
  primal_dual_pb *pdProb = setupProblem(prob);

#ifdef REP_MEM
  reportMem();
#endif

  // write the deterministic equivalent in mps format
  if (opt.writeMps()) {
    FILE *fout = fopen("smps.mps", "w");
    Write_MpsFile(fout, pdProb->AlgAug, pdProb->b, pdProb->c, pdProb->u, 0,
		  prob->colnames, prob->rownames);
    fclose(fout);
  }

  printf(" --------------- oopsSolve -----------------\n");

  // options for the complete problem
  hopdm_options->glopt->conv_tol = 1.e-4;

  hopdm_ret *ret;
  PrintOptions *Prt = NewHopdmPrt(PRINT_ITER);

  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // start the clock
  tt_start = oopstime();

  // solve the problem
  ret = hopdm(printout, pdProb, hopdm_options, Prt);

  // stop the clock
  tt_end = oopstime();

  // report statistics
  fprintf(printout, "Elapsed time: %.10g seconds.\n", tt_end - tt_start);

  if (opt.printSolution())
    getSolution(pdProb, prob);

  // clean up
  free(ret);

 TERMINATE:

  // clean up
  FreePDProblem(pdProb);
  freeSmpsReturn(prob);
  delete Prt;

  return 0;
}

/**
 *  Set up the Oops algebras and vectors and build the primal-dual problem.
 *
 *  @param Pb:
 *         Pointer to an SmpsReturn structure.
 *  @return Pointer to a primal_dual_pb structure containing the problem
 *          to be solved by Oops.
 *
 *  @note
 *  InitAlgebrasNew assumes that a callback function is set up for all
 *  SparseMatrix leaves of A and Q.
 */
primal_dual_pb* SmpsOops::setupProblem(SmpsReturn *Pb) {

  Algebra *A = Pb->AlgA;
  Algebra *Q = Pb->AlgQ;

  printf(" --------------- InitAlgebras --------------\n");
  Algebra *AlgAug = InitAlgebrasNew(A, Q);

  Vector *vb = NewVector(A->Trow, "vb");
  Vector *vc = NewVector(A->Tcol, "vc");
  Vector *vu = NewVector(A->Tcol, "vu");

  Vector *vx = NewVector(A->Tcol, "vx");
  Vector *vy = NewVector(A->Trow, "vy");
  Vector *vz = NewVector(A->Tcol, "vz");

  CopyDenseToVector(Pb->b, vb);
  CopyDenseToVector(Pb->c, vc);
  CopyDenseToVector(Pb->u, vu);

  return NewPDProblem(AlgAug, vb, vc, vu, vx, vy, vz);
}

/**
 *  Retrieve the solution.
 *
 *  @param Prob:
 *         Pointer to the problem solved by OOPS.
 *  @param Ret:
 *         Pointer to an SmpsReturn structure.
 *
 *  @bug
 *  The value of the slacks is always zero, as the slacks would need to be
 *  computed here.
 */
int SmpsOops::getSolution(primal_dual_pb *Prob, SmpsReturn *Ret) {

  DenseVector *x, *y, *z, *r;
  Tree *Trow = Ret->AlgA->Trow;
  Tree *Tcol = Ret->AlgA->Tcol;
  int nRows  = Trow->end - Trow->begin;
  int nCols  = Tcol->end - Tcol->begin;

  // allocate space for the solution vectors
  x = NewDenseVector(nCols, "dx");
  y = NewDenseVector(nRows, "dy");
  z = NewDenseVector(nCols, "dz");

  // this is just a placeholder, we cannot yet compute the slacks
  r = NewDenseVector(nRows, "slacks");

  // recover initial order on solution vectors
  SmpsVectorToDense(Prob->x, x, Ret, ORDER_COL);
  SmpsVectorToDense(Prob->y, y, Ret, ORDER_ROW);
  SmpsVectorToDense(Prob->z, z, Ret, ORDER_COL);

  // print the solution
  const NodeInfo *info = getNodeInfo();

#ifdef WITH_MPI
  if (IS_ROOT_PAR)
#endif
    printSolution(info, x->elts, y->elts, r->elts, z->elts);

  // clean up
  delete[] info->nRowsNode;
  delete[] info->nColsNode;
  delete info;

  FreeDenseVector(x);
  FreeDenseVector(y);
  FreeDenseVector(z);
  FreeDenseVector(r);

  return 0;
}

/** Return information on the nodes necessary for printing the solution */
const NodeInfo* SmpsOops::getNodeInfo() const {
  return smps.getNodeInfo();
}

/**
 *  Order nodes according to the level.
 *
 *  The nodes of the event tree are reordered according to the cutoff level:
 *  - breadth-first until given level (so that these nodes come first)
 *  - depth-first afterwards          (to give diagonal structure in rest)
 * Sets the following elements of tree:
 *  order[nd]     giving the original position (name of node) that
 *                      is at block position 'nd' in the big matrix
 *  block[nd]     giving the BLOCK that block 'nd' in big matrix
 *                      belongs to.
 *                      BLOCK = 0: root BLOCK (i.e. rankCor)
 *                      BLOCK = 1 -- TREE->nBlocks are diagonal parts
 */
int SmpsOops::orderNodes() {

  int nPeriods = smps.getPeriods();

  // check that the required cutoff level is feasible
  if (level <= 0 || level >= nPeriods) {
    printf("The cutoff level must be positive and smaller than "
           "the number of stages (%d).\n",  nPeriods);
    return 1;
  }

  // current node in the new list which is being followed
  int node = 0;

  // next node in the new list that is being built
  int nx_nd = 1;

  int nNodes = smps.getNodes();

  block    = new int[nNodes];
  order    = new int[nNodes];
  revorder = new int[nNodes];
  order[0] = 0;
  block[0] = 0;

  printf("Ordering %d nodes with %d levels in rnkcr.\n", nNodes, level);

  // make list of all nodes in rank corrector part (block 0)
  while (smps.getPeriod(order[node]) < level) {

    // add children of current node to list
    for (int i = 0; i < smps.getNChildren(order[node]); i++) {
      order[nx_nd] = smps.getFirstChild(order[node]) + i - 1;
      block[nx_nd] = 0;
      nx_nd++;
    }
    node++;
  }

  /* now all nodes up to period indicated by level are on the lists
      nx_nd points to the next free position in the list,
      node points to the first node of period = level in the list */

  // now follow depth first
  while (smps.getPeriod(order[node]) == level && node < nNodes) {

    // each node in this list is a seed for a diagonal block
    addChildrenToList(node, &nx_nd);
    node++;
  }

  // generate reverse order of nodes
  for (int i = 0; i < nNodes; i++)
    revorder[order[i]] = i;

#ifdef DEBUG_ORDER
  // report order of nodes
  printf("Found %d diagonal blocks.\n", nBlocks);
  printf("Reporting new order of nodes in big Matrix:\n");
  printf("   pos  node   per  block\n");
  for (int i = 0; i < nNodes; i++)
    printf("  %3d   %3d   %3d   %3d\n", i, order[i],
	   smps.getPeriod(order[i]), block[i]);
#endif

  // reset the period starts
  smps.setNodeStarts(order);

  return 0;
}

/** Add all children of node order[node] to the list starting at entry next */
void SmpsOops::addChildrenToList(const int node, int *next) {

  int i;
  const int ordNode = order[node];

  for (i = 0; i < smps.getNChildren(ordNode); i++) {
    if (smps.getPeriod(ordNode) == level)
      nBlocks++;

    order[*next] = smps.getFirstChild(ordNode) + i - 1;
    block[*next] = nBlocks;
    (*next)++;
    addChildrenToList((*next) - 1, next);
  }
}

/** Constructor */
OptionsOops::OptionsOops(const int argc, const char *argv[]) :
  Options(argc, argv),
  _useWarmstart(0),
  _cutoffLevel(1) {
}

/** Parse the command line options */
int OptionsOops::parse() {

  // add the specialised options
  Options::addOption("-w", "use a warmstart strategy",
		     &_useWarmstart);
  Options::addOption("-l", "cutoff level (for multistage programs)",
		     &_cutoffLevel, true);

  // parse the common options
  int rv = Options::parse();

  return rv;
}
