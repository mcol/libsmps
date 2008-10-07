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

#include <queue>
#include <assert.h>
#include "SmpsOops.h"
#include "oops/oopstime.h"
#include "oops/WriteMps.h"


double tt_start, tt_end;

static void
dfsMap(map<const Node*, Node*> &nMap, const Node *cNode, Node *rNode);


/** Constructor */
SmpsOops::SmpsOops(string smpsFile, const int lev) :
  smps(smpsFile),
  rTree(),
  pdPoint(NULL),
  wsPoint(NULL),
  cutoff(lev),
  nBlocks(0) {
}

/** Destructor */
SmpsOops::~SmpsOops() {

  delete pdPoint;
  delete wsPoint;
}

/** Read the smps files */
int SmpsOops::read() {

  int rv = smps.read(true);
  if (rv)
    return rv;

  // reorder the nodes according to the level
  rv = orderNodes(smps.getSmpsTree());

  return rv;
}

/**
 *  Solve the complete problem.
 *
 *  @param opt:
 *         Command line options.
 *  @param hopdmOpts:
 *         Options for the solver.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::solve(const OptionsOops &opt, HopdmOptions &hopdmOpts) {

  int rv = 0;

  // XXX This call is needed because nBlocks is in SmpsOops and we need it to
  // update it, since otherwise we may still use the one of the reduced tree
  if (rTree.getRootNode())
    orderNodes(smps.getSmpsTree());

  if (opt.writeMps()) {
    smps.setBuildNames();
  }

  SmpsReturn prob;

  // generate the problem
  rv = generateSmps(smps.getSmpsTree(), prob);
  if (rv) {
    printf("Failed to generate the deterministic equivalent.\n");
    return 1;
  }

  // setup the primal-dual problem
  PDProblem *pdProb = setupProblem(prob);

  // write the deterministic equivalent in mps format
  if (opt.writeMps()) {
    FILE *fout = fopen("smps.mps", "w");
    Write_MpsFile(fout, pdProb->AlgAug, pdProb->b, pdProb->c,
		  pdProb->u, pdProb->l, 0, prob.colnames, prob.rownames);
    fclose(fout);
  }

  printf(" --------------- solve ---------------------\n");

  hopdm_ret *ret = NULL;
  PrintOptions Prt(hopdmOpts.glopt->prt);

  // exit early if we don't have to solve the problem
  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // options for the complete problem
  hopdmOpts.glopt->conv_tol = 1.e-4;

  // use the warmstart point if available
  if (wsPoint)
    hopdmOpts.use_start_point = 1;

  // start the clock
  tt_start = oopstime();

  // solve the problem
  ret = hopdm(printout, pdProb, &hopdmOpts, &Prt);
  if (ret->ifail) {
    rv = ret->ifail;
    goto TERMINATE;
  }

  // stop the clock
  tt_end = oopstime();

  // report statistics
  fprintf(printout, "Elapsed time: %.10g seconds.\n", tt_end - tt_start);

  if (opt.printSolution())
    getSolution(pdProb, prob);

 TERMINATE:

  // clean up
  delete ret;
  FreePDProblem(pdProb);

  return rv;
}

/**
 *  Solve a reduced problem.
 *
 *  @param opt:
 *         Command line options.
 *  @param hopdmOpts:
 *         Options for the solver.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::solveReduced(const OptionsOops &opt,
			   HopdmOptions &hopdmOpts) {

  int rv = 0;

  SmpsReturn prob;

  // generate a reduced problem
  rv = generateSmps(rTree, prob);
  if (rv) {
    printf("Failed to generate the deterministic equivalent.\n");
    return 1;
  }

  // setup the primal-dual problem
  PDProblem *pdProb = setupProblem(prob);

  // write the deterministic equivalent in mps format
  if (opt.writeMps()) {
    FILE *fout = fopen("smps-red.mps", "w");
    Write_MpsFile(fout, pdProb->AlgAug, pdProb->b, pdProb->c,
		  pdProb->u, pdProb->l, 0, prob.colnames, prob.rownames);
    fclose(fout);
  }

  printf(" --------------- solveReduced --------------\n");

  hopdm_ret *ret = NULL;
  PrintOptions Prt(hopdmOpts.glopt->prt);

  // exit early if we don't have to solve the problem
  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    goto TERMINATE;
  }

  // options for the reduced problem
  hopdmOpts.glopt->conv_tol = 5.e-1;

  // solve the problem
  ret = hopdm(printout, pdProb, &hopdmOpts, &Prt);
  if (ret->ifail) {
    rv = ret->ifail;
    goto TERMINATE;
  }

  // extract and store the solution
  storeSolution(pdProb, prob);

  // generate a warmstart point for the complete problem
  setupWarmStart(prob);

  if (pdPoint) {
    FreeVector(pdPoint->x);
    FreeVector(pdPoint->y);
    FreeVector(pdPoint->z);
    if (smps.hasUpperBounds()) {
      FreeVector(pdPoint->s);
      FreeVector(pdPoint->w);
    }
  }

 TERMINATE:

  // clean up
  delete ret;
  FreePDProblem(pdProb);

  return rv;
}

/**
 *  Store the solution from the reduced problem.
 *
 *  @param pdProb:
 *         The reduced primal-dual problem.
 *  @param Ret:
 *         The SmpsReturn structure of the reduced problem.
 *  @return 1 If something goes wrong, 0 otherwise.
 */
int SmpsOops::storeSolution(const PDProblem *pdProb, const SmpsReturn &Ret) {

  // allocate space for the new vectors
  Algebra *A = Ret.AlgA;
  Vector *x = NewVector(A->Tcol, "vx");
  Vector *y = NewVector(A->Trow, "vy");
  Vector *z = NewVector(A->Tcol, "vz");
  Vector *s = NULL, *w = NULL;
  if (smps.hasUpperBounds()) {
    s = NewVector(A->Tcol, "vs");
    w = NewVector(A->Tcol, "vw");
  }

  // copy the solution vectors
  CopyVector(pdProb->x, x);
  CopyVector(pdProb->y, y);
  CopyVector(pdProb->z, z);
  if (smps.hasUpperBounds()) {
    CopyVector(pdProb->s, s);
    CopyVector(pdProb->w, w);
  }

  // store the solution
  pdPoint = NewPDPoint(x, y, z, s, w, NULL);

  return 0;
}

/**
 *  Create the reduced tree in a recursive manner.
 *
 *  @param cNode:
 *         Node in the complete tree.
 *  @param rParent:
 *         Node in the reduced tree to which add children.
 *  @param nWanted:
 *         Number of desired scenarios in the reduced tree.
 */
void SmpsOops::reduceScenarios(const Node *cNode, Node *rParent,
			       const int nWanted) {

  int i, nChildren = cNode->nChildren();
  int chdIdx = 0, nAdded = 0;
  int each = 0, rest = 0, step = 1;
  Node *child;

  // for each stage cache the reduced tree node from which the
  // complete tree nodes should be initialised
  static Node* redChd[MAX_PERIODS];

#ifdef DEBUG_RTREE
  printf("ReduceScenarios: node %d (%d)\n", cNode->name(), nWanted);
#endif

  if (nChildren > 0) {
    // common number of scenarios to be chosen from each child
    each = nWanted / nChildren;

    // number of nodes that have an extra scenario
    rest = nWanted % nChildren;
  }

  if (nWanted > 0)
    step = nChildren / nWanted;

  // copy the needed number of children of the complete node
  // we copy only the nodes that are at position chdIdx
  for (i = 0; i < nChildren; ++i) {

    Node *ttt = cNode->getChild(i);

    // this node is not going in the reduced tree
    if (i != chdIdx || nWanted == 0) {
      assert(redChd[ttt->level()] != NULL);

      // just map the node and its children
      nMap[ttt] = redChd[ttt->level()];
      reduceScenarios(ttt, child, each);
      continue;
    }

    // copy the node in the reduced tree
    child = new Node(ttt->name());
    child->copy(ttt);
    rParent->addChild(child);
    nMap[ttt] = child;

    // store the child from which to initialise the nodes that
    // will not be put in the reduced tree
    redChd[ttt->level()] = child;

    // index of the next child to be put in the reduced tree
    chdIdx += step;

    // stop adding nodes to the reduced tree
    if (++nAdded >= nWanted)
      chdIdx = -1;

#ifdef DEBUG_RTREE
    printf("Node %d: chosen for the reduced tree.\n", ttt->name());
#endif

    // enter the recursion
    if (rest > 0) {
      reduceScenarios(ttt, child, each + 1);
      --rest;
    }
    else
      reduceScenarios(ttt, child, each);
  }
}

/**
 *  Generate a reduced tree by choosing a subset of scenarios.
 *
 *  @param nScenarios:
 *         Number of scenarios to appear in the reduced tree.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::reduceTree(const int nScenarios) {

  printf(" --------------- reduceTree ----------------\n");

  // check that the number of scenarios is sensible
  if (nScenarios < 1)
    return 1;

  const Node *cNode = smps.getRootNode();

  // allocate the root node for the reduced tree
  rTree.setRootNode(new Node(100 + cNode->name()));
  Node *rNode = rTree.getRootNode();

  // copy the root node
  rNode->copy(cNode);
  nMap[cNode] = rNode;

  // build up the reduced tree by selecting some scenarios
  reduceScenarios(cNode, rNode, nScenarios);

  // set the start rows and columns for each node
  int rv = smps.setNodeStarts(rTree);
  if (rv)
    return 1;

  // order the nodes and set the next links
  orderNodes(rTree);

  // recompute the probabilities in the reduced tree
  adjustProbabilities();

#ifdef DEBUG_RTREE
  rTree.print();
#endif

  return 0;
}

/**
 *  Generate a reduced tree with aggregation.
 *
 *  @param nAggr:
 *         Number of stages to be aggregated (starting from the last)
 *  @return 1 If something goes wrong, 0 otherwise.
 */
int SmpsOops::aggregateStages(const int nAggr) {

  // check that the number of stages to aggregate is sensible
  if (nAggr < 1)
    return 1;

  printf(" --------------- aggregateStages -----------\n");

  // aggregate the last stages
  const int last = smps.getPeriods() - nAggr + 1;
  if (last < 2) {
    printf("No aggregation possible.\n");
    return 1;
  }

  printf("Aggregating stages %d to %d.\n", last, smps.getPeriods());

  int nodeName = 100;
  const Node *cNode = smps.getRootNode();
  Node *rNode;

  // queue of nodes to be processed
  queue<const Node*> cNodes;
  cNodes.push(cNode);

  // copy the tree up to the aggregation point
  while (!cNodes.empty()) {

    // take the first element in the queue
    cNode = cNodes.front();
    cNodes.pop();

    // these nodes have to be copied
    if (cNode->level() < last) {

      // create a node in the reduced tree
      rNode = new Node(++nodeName);
      rNode->copy(cNode);
      rNode->setProb(cNode->probNode());
      if (cNode->parent())
        nMap[cNode->parent()]->addChild(rNode);

      // map the complete node to the reduced node
      nMap[cNode] = rNode;

      // these nodes span from the current stage to the last
      if (rNode->level() == last - 1) {
	rNode->setLevel(nAggr);

	// map all children to this reduced node
	dfsMap(nMap, cNode, rNode);
      }
    }

    // put the children in the queue
    for (int i = 0; i < cNode->nChildren(); ++i)
      cNodes.push(cNode->getChild(i));
  }

  // find the root node
  while (rNode->level() > 0)
    rNode = rNode->parent();

  // set the root of the reduced tree
  rTree.setRootNode(rNode);

  // order the nodes and set the next links
  orderNodes(rTree);

#ifdef DEBUG_RTREE
  rTree.print();
#endif

  return 0;
}

void dfsMap(map<const Node*, Node*> &nMap, const Node *cNode, Node *rNode) {

  const Node *child;
  for (int i = 0; i < cNode->nChildren(); ++i) {
    child = cNode->getChild(i);
    nMap[child] = rNode;
    dfsMap(nMap, child, rNode);
  }
}

/**
 *  Set up the Oops algebras and vectors and build the primal-dual problem.
 *
 *  @param Pb:
 *         The SmpsReturn structure of the problem generated.
 *  @return Pointer to the problem to be solved by Oops.
 *
 *  @note
 *  InitAlgebrasNew assumes that a callback function is set up for all
 *  SparseMatrix leaves of A and Q.
 */
PDProblem* SmpsOops::setupProblem(SmpsReturn &Pb) {

  Algebra *A = Pb.AlgA;
  Algebra *Q = Pb.AlgQ;

  Algebra *AlgAug = InitAlgebrasNew(A, Q);

  Vector *vb = NewVector(A->Trow, "vb");
  Vector *vc = NewVector(A->Tcol, "vc");
  Vector *vu = NewVector(A->Tcol, "vu");
  Vector *vl = NewVector(A->Tcol, "vl");

  Vector *vx = NewVector(A->Tcol, "vx");
  Vector *vy = NewVector(A->Trow, "vy");
  Vector *vz = NewVector(A->Tcol, "vz");
  Vector *vs = NULL, *vw = NULL;
  if (smps.hasUpperBounds()) {
    vs = NewVector(A->Tcol, "vs");
    vw = NewVector(A->Tcol, "vw");
  }

  // use the warmstart point if available
  if (wsPoint) {
    SmpsDenseToVector(wsPoint->x, vx, Pb, ORDER_COL);
    SmpsDenseToVector(wsPoint->y, vy, Pb, ORDER_ROW);
    SmpsDenseToVector(wsPoint->z, vz, Pb, ORDER_COL);
    if (smps.hasUpperBounds()) {
      SmpsDenseToVector(wsPoint->s, vs, Pb, ORDER_COL);
      SmpsDenseToVector(wsPoint->w, vw, Pb, ORDER_COL);
    }
  }

  CopyDenseToVector(Pb.b, vb);
  CopyDenseToVector(Pb.c, vc);
  CopyDenseToVector(Pb.u, vu);
  CopyDenseToVector(Pb.l, vl);

  // create the primal dual problem
  PDProblem *Prob = NewPDProblem(AlgAug, vb, vc, vu, vx, vy, vz);
  if (smps.hasUpperBounds()) {
    Prob->s = vs;
    Prob->w = vw;
  }
  Prob->l= vl;

  return Prob;
}

/** Recompute the probabilities in the reduced tree */
void SmpsOops::adjustProbabilities() {

  const Node *cNode = smps.getRootNode();

  do {

    // find the corresponding node in the reduced tree
    Node *rNode = nMap[cNode];

    // update the probabilities
    rNode->setProb(rNode->probNode() + cNode->probNode());

  } while (cNode = cNode->next());
}

/**
 *  Set up a warmstart point from a reduced-tree solution.
 *
 *  @param Ret:
 *         The SmpsReturn structure of the reduced problem.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::setupWarmStart(const SmpsReturn &Ret) {

  // dense vectors for the reduced solution
  DenseVector *xred, *zred, *yred, *sred = NULL, *wred = NULL;

  // dense vectors for the complete solution
  DenseVector *xnew, *znew, *ynew, *snew = NULL, *wnew = NULL;

  printf(" --------------- setupWarmStart ------------\n");

  if (!pdPoint) {
    printf("Failed to set up a warmstart point.\n");
    return 1;
  }

  // dimensions of the complete and the reduced deterministic equivalents
  int nRows = smps.getTotRows(), rRows = rTree.getTotRows();
  int nCols = smps.getTotCols(), rCols = rTree.getTotCols();

  printf("Reduced matrix:  %dx%d\n", rRows, rCols);
  printf("Complete matrix: %dx%d\n", nRows, nCols);

#ifdef DEBUG_RTREE
  // print the mapping between complete and reduced tree
  map<const Node*, Node*>::iterator it;
  printf("Mapping: cNode => rNode:\n");
  for (it = nMap.begin(); it != nMap.end(); it++)
    printf("%3d => %3d\n", (*it).first->name(), (*it).second->name());
#endif

  // allocate space for the vectors in the reduced iterate
  xred = NewDenseVector(rCols, "xred");
  yred = NewDenseVector(rRows, "yred");
  zred = NewDenseVector(rCols, "zred");
  if (smps.hasUpperBounds()) {
    sred = NewDenseVector(rCols, "sred");
    wred = NewDenseVector(rCols, "wred");
  }

  // recover the initial ordering of the solution vectors
  VectorToSmpsDense(pdPoint->x, xred, Ret, ORDER_COL);
  VectorToSmpsDense(pdPoint->y, yred, Ret, ORDER_ROW);
  VectorToSmpsDense(pdPoint->z, zred, Ret, ORDER_COL);
  if (smps.hasUpperBounds()) {
    VectorToSmpsDense(pdPoint->s, sred, Ret, ORDER_COL);
    VectorToSmpsDense(pdPoint->w, wred, Ret, ORDER_COL);
  }

  // allocate space for the vectors in the complete iterate
  xnew = NewDenseVector(nCols, "xnew");
  ynew = NewDenseVector(nRows, "ynew");
  znew = NewDenseVector(nCols, "znew");
  if (smps.hasUpperBounds()) {
    snew = NewDenseVector(nCols, "snew");
    wnew = NewDenseVector(nCols, "wnew");
  }

  const Node *cNode = smps.getRootNode(), *rNode;
  int cIndex, rIndex, nElems;
  double crProb;

  // go through the nodes in the complete tree
  do {

    rNode  = nMap[cNode];
    cIndex = cNode->firstCol();
    rIndex = rNode->firstCol();
    nElems = cNode->nCols();
    crProb = cNode->probNode() / rNode->probNode();

    // adjust rIndex taking into account aggregation
    if (cNode->level() != rNode->level()) {
      const Node *ttt = cNode;
      while (ttt->level() != rNode->level()) {
	ttt = ttt->parent();
	rIndex += ttt->nCols();
      }
    }

#ifdef DEBUG_WARMSTART
    printf("Working on %d -> rNode: %d,  %d  %d\n",
	   cNode->name(), rNode->name(), cIndex, rIndex);
#endif

    // copy the x and z iterates
    for (int i = 0; i < nElems; i++) {
      xnew->elts[cIndex + i] = xred->elts[rIndex + i];
      znew->elts[cIndex + i] = zred->elts[rIndex + i] * crProb;
      if (smps.hasUpperBounds()) {
	snew->elts[cIndex + i] = sred->elts[rIndex + i];
	wnew->elts[cIndex + i] = wred->elts[rIndex + i] * crProb;
      }
    }

    cIndex = cNode->firstRow();
    rIndex = rNode->firstRow();
    nElems = cNode->nRows();

    // adjust rIndex taking into account aggregation
    if (cNode->level() != rNode->level()) {
      const Node *ttt = cNode;
      while (ttt->level() != rNode->level()) {
	ttt = ttt->parent();
	rIndex += ttt->nRows();
      }
    }

    // copy the y iterate
    for (int i = 0; i < nElems; i++) {
      ynew->elts[cIndex + i] = yred->elts[rIndex + i] * crProb;
    }

  } while (cNode = cNode->next());

  // set the warmstart point
  wsPoint = new WSPoint(xnew, ynew, znew, snew, wnew);

  // clean up
  FreeDenseVector(xred);
  FreeDenseVector(yred);
  FreeDenseVector(zred);
  if (smps.hasUpperBounds()) {
    FreeDenseVector(sred);
    FreeDenseVector(wred);
  }

  return 0;
}

/**
 *  Retrieve the solution.
 *
 *  @param Prob:
 *         Pointer to the problem solved by OOPS.
 *  @param Ret:
 *         The SmpsReturn structure of the problem.
 *
 *  @bug
 *  The value of the slacks is always zero, as the slacks would need to be
 *  computed here.
 */
int SmpsOops::getSolution(PDProblem *Prob, SmpsReturn &Ret) {

  DenseVector *x, *y, *z, *r;
  Tree *Trow = Ret.AlgA->Trow;
  Tree *Tcol = Ret.AlgA->Tcol;
  int nRows  = Trow->end - Trow->begin;
  int nCols  = Tcol->end - Tcol->begin;

  // allocate space for the solution vectors
  x = NewDenseVector(nCols, "dx");
  y = NewDenseVector(nRows, "dy");
  z = NewDenseVector(nCols, "dz");

  // this is just a placeholder, we cannot yet compute the slacks
  r = NewDenseVector(nRows, "slacks");

  // recover initial order on solution vectors
  VectorToSmpsDense(Prob->x, x, Ret, ORDER_COL);
  VectorToSmpsDense(Prob->y, y, Ret, ORDER_ROW);
  VectorToSmpsDense(Prob->z, z, Ret, ORDER_COL);

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
int SmpsOops::orderNodes(SmpsTree &Tree) {

  int nPeriods = smps.getPeriods();
  Node *node = Tree.getRootNode();

  // leave immediately if there is no root node
  if (!node)
    return 1;

  // check that the required cutoff level is feasible
  if (cutoff <= 0 || cutoff >= nPeriods) {
    printf("The cutoff level must be positive and smaller than "
           "the number of stages (%d).\n", nPeriods);
    return 1;
  }

  // reset the cutoff level to the maximum supported
  if (cutoff > MAX_CUTOFF)
    cutoff = MAX_CUTOFF;

  // reset the number of blocks, because we may call orderNodes multiple
  // times and on different trees
  nBlocks = 0;

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

    // we should never have an empty queue
    assert(!qNodes.empty());

    // take the first node in the queue
    node = qNodes.front();

  } while (node->level() < cutoff);

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
  node = Tree.getRootNode();

  do {

    qOrder.pop();
    node->setNext(qOrder.empty() ? NULL : qOrder.front());

  } while (node = node->next());

#ifdef DEBUG_ORDER
  // report order of nodes
  printf("Found %d diagonal blocks.\n", nBlocks);
  printf("Reporting new order of nodes in big Matrix:\n");
  printf("  node   per  block\n");
  node = Tree.getRootNode();
  do {
    printf("  %3d   %3d   %3d\n",
	   node->name(), node->level() + 1, node->block());
  } while (node = node->next());
#endif

  // reset the period starts
  smps.setNodeStarts(Tree);

  // dimensions of the deterministic equivalent
  printf("The deterministic equivalent matrix is %dx%d, %d nonzeros.\n",
	 Tree.getTotRows(), Tree.getTotCols(), smps.countNonzeros(Tree));

  return 0;
}

/** Perform a recursive depth-first ordering of the node and its children */
void SmpsOops::dfsNode(queue<Node*> &qOrder, Node *node) {

  // put the node in the reordered queue
  qOrder.push(node);

  // count the number of diagonal blocks
  if (node->level() == cutoff)
    ++nBlocks;

  node->setBlock(nBlocks);

  // traverse the children in depth-first order
  for (int i = 0; i < node->nChildren(); i++) {
    dfsNode(qOrder, node->getChild(i));
  }
}

/** Constructor */
SmpsReturn::SmpsReturn() :
  AlgA(NULL),
  AlgQ(NULL),
  b(NULL),
  c(NULL),
  l(NULL),
  u(NULL),
  rownames(NULL),
  colnames(NULL),
  rootNode(NULL),
  is_col_diag(NULL),
  nRowsRnkc(0),
  nColsRnkc(0),
  nColsDiag(0) {
}

/** Destructor */
SmpsReturn::~SmpsReturn() {

  if (AlgA) {
    FreeTree(AlgA->Trow);
    FreeTree(AlgA->Tcol);
    FreeAlgebraAlg(AlgA);
  }

  if (AlgQ) {
    FreeTree(AlgQ->Trow);
    FreeTree(AlgQ->Tcol);
    FreeAlgebraAlg(AlgQ);
  }

  if (rownames) {
    const int ttm = b->dim;
    for (int i = 0; i < ttm; ++i)
      delete[] rownames[i];
    delete[] rownames;
  }

  if (colnames) {
    const int ttn = c->dim;
    for (int i = 0; i < ttn; ++i)
      delete[] colnames[i];
    delete[] colnames;
  }

  FreeDenseVector(b);
  FreeDenseVector(c);
  FreeDenseVector(l);
  FreeDenseVector(u);

  delete[] is_col_diag;
}

/** Constructor */
WSPoint::WSPoint(DenseVector *vx, DenseVector *vy, DenseVector *vz,
		 DenseVector *vs, DenseVector *vw) :
  x(vx),
  y(vy),
  z(vz),
  s(vs),
  w(vw) {
}

/** Destructor */
WSPoint::~WSPoint() {

  FreeDenseVector(x);
  FreeDenseVector(y);
  FreeDenseVector(z);
  if (s)
    FreeDenseVector(s);
  if (w)
    FreeDenseVector(w);
}

/** Constructor */
OptionsOops::OptionsOops(const int argc, const char *argv[]) :
  Options(argc, argv),
  _useReduction(0),
  _useAggregation(0),
  _cutoffLevel(1) {
}

/** Parse the command line options */
int OptionsOops::parse() {

  // add the specialised options
  Options::addOption("-w", "use a warmstart strategy with scenario reduction",
		     &_useReduction, true);
  Options::addOption("-a", "use a warmstart strategy with stage aggregation",
		     &_useAggregation, true);
  Options::addOption("-l", "cutoff level (for multistage programs)",
		     &_cutoffLevel, true);

  // parse the common options
  int rv = Options::parse();

  return rv;
}
