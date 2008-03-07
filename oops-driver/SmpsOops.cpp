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
#include <queue>
#include "SmpsOops.h"
#include "oops/oopstime.h"
#include "oops/WriteMps.h"


double tt_start, tt_end;


/** Constructor */
SmpsOops::SmpsOops(string smpsFile, const int lev) :
  smps(smpsFile),
  level(lev),
  nBlocks(0) {
}

/** Destructor */
SmpsOops::~SmpsOops() {

}

/** Read the smps files */
int SmpsOops::read() {

  int rv = smps.read();
  if (rv)
    return rv;

  // take care of slacks
  smps.modifyCore();

  // reorder the nodes according to the level
  rv = orderNodes(smps.getRootNode());

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
  PDProblem *pdProb = setupProblem(prob);

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
 *  @return Pointer to the problem to be solved by Oops.
 *
 *  @note
 *  InitAlgebrasNew assumes that a callback function is set up for all
 *  SparseMatrix leaves of A and Q.
 */
PDProblem* SmpsOops::setupProblem(SmpsReturn *Pb) {

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
int SmpsOops::getSolution(PDProblem *Prob, SmpsReturn *Ret) {

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
#ifdef WITH_MPI
  if (IS_ROOT_PAR)
#endif
    printSolution(smps.getRootNode(), x->elts, y->elts, r->elts, z->elts);

  // clean up
  FreeDenseVector(x);
  FreeDenseVector(y);
  FreeDenseVector(z);
  FreeDenseVector(r);

  return 0;
}

/**
 *  Order nodes according to the level.
 *
 *  The nodes of the event tree are reordered according to the cutoff level:
 *  - breadth-first until given level (so that these nodes come first)
 *  - depth-first afterwards          (to give diagonal structure in rest)
 *
 *  @param node:
 *         Root node of the tree to be reordered.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::orderNodes(Node *rootNode) {

  int nPeriods = smps.getPeriods();
  Node *node = rootNode;

  // leave immediately if there is no root node
  if (!node)
    return 1;

  // check that the required cutoff level is feasible
  if (level <= 0 || level >= nPeriods) {
    printf("The cutoff level must be positive and smaller than "
           "the number of stages (%d).\n", nPeriods);
    return 1;
  }

  // queue of nodes to be followed
  queue<Node*> qNodes;

  // queue of nodes in the new order
  queue<Node*> qOrder;

  // start from the root
  qNodes.push(node);

  // put the nodes up to the cutoff level in breadth-first order
  do {

    // take the first element from the queue
    qNodes.pop();

    // put the node in the reordered list
    qOrder.push(node);

    // add its children to the queue
    for (int i = 0; i < node->nChildren(); ++i) {
      qNodes.push(node->getChild(i));
    }

    // take the first node in the queue
    node = qNodes.front();

  } while (node->level() < level);

  // now proceed in depth-first order
  while (!qNodes.empty()) {

    // take the first element in the queue
    node = qNodes.front();
    qNodes.pop();

    // call a recursive function
    dfsNode(qOrder, node);
  }

  // update the next links

  // start from the root
  node = rootNode;

  do {

    qOrder.pop();
    node->setNext(qOrder.empty() ? NULL : qOrder.front());

  } while (node = node->next());

#ifdef DEBUG_ORDER
  // report order of nodes
  printf("Found %d diagonal blocks.\n", nBlocks);
  printf("Reporting new order of nodes in big Matrix:\n");
  printf("  node   per  block\n");
  node = rootNode;
  do {
    printf("  %3d   %3d   %3d\n",
	   node->name(), node->level() + 1, node->block());
  } while (node = node->next());
#endif

  // reset the period starts
  smps.setNodeStarts(rootNode);

  return 0;
}

/** Perform a recursive depth-first ordering of the node and its children */
void SmpsOops::dfsNode(queue<Node*> &qOrder, Node *node) {

  // put the node in the reordered queue
  qOrder.push(node);

  // count the number of diagonal blocks
  if (node->level() == level)
    ++nBlocks;

  node->setBlock(nBlocks);

  // traverse the children in depth-first order
  for (int i = 0; i < node->nChildren(); i++) {
    dfsNode(qOrder, node->getChild(i));
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
