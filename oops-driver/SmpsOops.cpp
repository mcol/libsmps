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
#include <cstring>
#include <assert.h>
#include "SmpsOops.h"


static void
dfsMap(map<const Node*, Node*> &nMap, const Node *cNode, Node *rNode);

static void
dfsNode(queue<Node*> &qOrder, Node *node, int &nBlocks, const int cutoff);


/** Constructor */
SmpsOops::SmpsOops(string smpsFile, const int lev) :
  fsContr(NULL),
  smps(smpsFile),
  rTree(),
  wsPoint(NULL),
  wsReady(false),
  cutoff(lev) {
}

/** Destructor */
SmpsOops::~SmpsOops() {

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
 *  @return A nonzero value if something goes wrong; 0 otherwise.
 */
int SmpsOops::solve(const OptionsOops &opt, HopdmOptions &hopdmOpts) {

  if (opt.writeMps()) {
    smps.setBuildNames();
  }

  // use the warmstart point if available
  if (wsReady)
    hopdmOpts.use_start_point = 1;

  printf(" --------------- solve ---------------------\n");

  // pass the problem to the solver
  int rv = solver(smps.getSmpsTree(), wsPoint, opt, hopdmOpts);

  return rv;
}

/**
 *  Solve a reduced problem.
 *
 *  @param opt:
 *         Command line options.
 *  @param hopdmOpts:
 *         Options for the solver.
 *  @return A nonzero value if something goes wrong; 0 otherwise.
 */
int SmpsOops::solveReduced(const OptionsOops &opt,
			   HopdmOptions &hopdmOpts) {

  // store the original convergence tolerance
  const double origTol = hopdmOpts.glopt->conv_tol;

  // options for the reduced problem
  hopdmOpts.glopt->conv_tol = 5.e-1;

  printf(" --------------- solveReduced --------------\n");

  // pass the problem to the solver
  int rv = solver(rTree, NULL, opt, hopdmOpts);
  if (rv)
    return rv;

  // the warmstart point is ready
  wsReady = true;

  // restore the tolerance
  hopdmOpts.glopt->conv_tol = origTol;

  return 0;
}

/**
 *  Solve a series of problems by decomposition.
 *
 *  It first generates and solves a problem based on a reduced tree in order
 *  to obtain an approximation to the first-stage decisions, which is then
 *  substituted out in the decomposed subproblems. If requested, these can be
 *  warmstarted from the solution obtained from the reduced tree.
 *
 *  @param opt:
 *         Command line options.
 *  @param hopdmOpts:
 *         Options for the solver.
 *  @return A nonzero value if something goes wrong; 0 otherwise.
 */
int SmpsOops::solveDecomposed(const OptionsOops &opt,
                              HopdmOptions &hopdmOpts) {

  printf(" --------------- solveDecomposed -----------\n");

  // check that the problem has at least two stages
  if (smps.getPeriods() < 2) {
    printf("No decomposition possible.\n");
    return 1;
  }

  // generate a reduced problem containing as many scenarios as there are
  // branches in the first stage unless a different value has been specified
  // from the command line
  int rv = reduceTree(opt.useReduction() ?
                      opt.useReduction() : smps.getRootNode()->nChildren());
  if (rv)
    return rv;

  rv = solveReduced(opt, hopdmOpts);
  if (rv)
    return rv;

  // reset the reduced tree to be empty
  rTree.reset();

  printf(" ---- Calling firstStageBlock ----\n");
  wsReady = false;

  // store the original convergence tolerance
  const double origTol = hopdmOpts.glopt->conv_tol;

  // options for the decomposed problem
  hopdmOpts.glopt->conv_tol = 5.e-2;

  // solve a subproblem rooted for each second-stage node
  const Node *root = smps.getRootNode();
  const int nSubProblems = root->nChildren();
  const bool bounds = smps.hasUpperBounds();

  // dense vectors for the subproblem warmstart
  WSPoint *wsDecomp = NULL;

  // compute the correction for the first stage variables
  fsContr = firstStageContribution();
  if (!fsContr) {
    rv = 1;
    goto TERMINATE;
  }

  // set up and solve the subproblems
  for (int chd = 0; chd < nSubProblems; ++chd) {

    printf(" --------------- subproblem %2d of %2d -------\n",
           chd + 1, nSubProblems);

    const Node *cNode = root->getChild(chd);

    // generate the subtree corresponding to the current child
    createSubtree(cNode, 1000 * (chd + 1));

    if (opt.warmstartDecomposition()) {

      // force an update of information held by the reduced tree so that the
      // dimension of the corresponding determinitic equivalent are correct
      orderNodes(rTree);

      // warmstart point for the subproblem
      wsDecomp = new WSPoint(rTree.getTotRows(), rTree.getTotCols(), bounds);

      // set up the warmstart point from the complete subtree rooted at cNode
      setupWSPoint(wsDecomp, wsPoint, cNode, true);
      hopdmOpts.use_start_point = 1;
      wsReady = true;
    }

    // pass the problem to the solver
    rv = solver(rTree, wsDecomp, opt, hopdmOpts);

    // clean up
    delete wsDecomp;
    rTree.reset();

    if (rv)
      goto TERMINATE;
  }

  // the warmstart point is now completely defined
  wsReady = true;

 TERMINATE:

  // restore the tolerance
  hopdmOpts.glopt->conv_tol = origTol;

  delete[] fsContr;
  fsContr = NULL;

  return rv;
}

/**
 *  Generate and solve the deterministic equivalent.
 *
 *  @param tree:
 *         The tree for the deterministic equivalent that is being built.
 *  @param wsIterate:
 *         The WSPoint to be used in case of warmstart.
 *  @param opt:
 *         Command line options.
 *  @param hopdmOpts:
 *         Options for the solver.
 *  @return A nonzero value if something goes wrong; 0 otherwise.
 */
int SmpsOops::solver(SmpsTree &tree, const WSPoint *wsIterate,
                     const OptionsOops &opt, HopdmOptions &hopdmOpts) {

  SmpsReturn prob;
  const bool reduced = &tree != &smps.getSmpsTree();

  if (reduced) {
    // order the nodes and set the next links
    orderNodes(tree);

#ifdef DEBUG_RTREE
    tree.print();
#endif
  }

  // generate the deterministic equivalent problem
  int rv = generateSmps(tree, prob);
  if (rv) {
    printf("Failed to generate the deterministic equivalent.\n");
    return rv;
  }

  // setup the primal-dual problem
  PDProblem pdProb = setupProblem(prob, wsIterate);

  // write the deterministic equivalent in mps format
  if (opt.writeMps()) {
    FILE *fout = fopen(reduced ? "smps-red.mps" : "smps.mps", "w");
    Write_MpsFile(fout, pdProb.AlgAug, pdProb.b, pdProb.c,
		  pdProb.u, pdProb.l, 0, prob.colnames, prob.rownames);
    fclose(fout);
  }

  // exit early if we don't have to solve the problem
  if (opt.dontSolve()) {
    printf("Problem not solved by request.\n");
    return 0;
  }

  PrintOptions Prt(reduced ? PRINT_NONE : PRINT_ITER);

  // solve the problem
  hopdm_ret *ret = hopdm(printout, &pdProb, &hopdmOpts, &Prt);
  rv = ret->ifail;
  free(ret);
  if (rv)
    return rv;

  // use the solution to a reduced problem to generate a warmstart point
  // for the complete problem
  if (reduced)
    setupWarmStart(pdProb, prob);

  // print the solution only for the complete problem
  if (opt.printSolution() && !reduced)
    getSolution(pdProb, prob);

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
  Node *cChild = NULL, *rChild = NULL;

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
  if (step == 0)
    step = 1;

  // copy the needed number of children of the complete node
  // we copy only the nodes that are at position chdIdx
  for (i = 0; i < nChildren; ++i) {

    cChild = cNode->getChild(i);

    // this node is not going in the reduced tree
    if (i != chdIdx || nWanted == 0) {
      assert(redChd[cChild->level()] != NULL);

      // just map the node and its children
      nMap[cChild] = redChd[cChild->level()];
      reduceScenarios(cChild, rChild, each);
      continue;
    }

    // copy the node in the reduced tree
    rChild = new Node(cChild->name());
    rChild->copy(cChild);
    rParent->addChild(rChild);
    nMap[cChild] = rChild;

    // store the child from which to initialise the nodes that
    // will not be put in the reduced tree
    redChd[cChild->level()] = rChild;

    // index of the next child to be put in the reduced tree
    chdIdx += step;

    // stop adding nodes to the reduced tree
    if (++nAdded >= nWanted)
      chdIdx = -1;

#ifdef DEBUG_RTREE
    printf("Node %d: chosen for the reduced tree.\n", cChild->name());
#endif

    // enter the recursion
    if (rest > 0) {
      reduceScenarios(cChild, rChild, each + 1);
      --rest;
    }
    else
      reduceScenarios(cChild, rChild, each);
  }
}

/**
 *  Create a reduced tree from a clustering information file.
 *
 *  Creates the reduced tree by using a clustering algorithm. The cluster
 *  information will be read in from an external file.
 *
 *  @note
 *  The first line of the clustering file is always skipped: the scenario
 *  mappings must start from the second line.
 *
 *  @note
 *  At the moment it works only for two stage problems.
 *
 *  @param clusteringFile:
 *         Name of the clustering file for scenario reduction.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::reduceScenariosCluster(const char *clusteringFile) {

  if (smps.getPeriods() > 2) {
    printf("Scenario clustering can be used only for two-stage problems.\n");
    return 1;
  }

  FILE *fin = fopen(clusteringFile, "r");
  if (!fin) {
    printf("Cannot find the clustering information file '%s'\n",
           clusteringFile);
    return 1;
  }

  // count the number of lines in the file
  const int nScenarios = smps.getMaxScens();
  int *clusters = new int[nScenarios], scen, src, nlines = 0;
  char buffer[100];

  // read the clustering file
  fgets(buffer, 100, fin);
  fgets(buffer, 100, fin);
  while (!feof(fin)) {
    nlines++;
    fgets(buffer, 100, fin);
  }

  // the clustering file must have as many lines as there are scenarios
  printf("Clustering file: %s (%d lines).\n", clusteringFile, nlines);
  if (nlines != nScenarios) {
    printf("Expected %d lines.\n", nScenarios);
    fclose(fin);
    delete[] clusters;
    return 1;
  }

  const Node *cNode = smps.getRootNode();
  Node *rParent = nMap[cNode];
  Node *child = NULL, *ttt;
  map<int, Node *> keepnodes;

  // now read the file again
  rewind(fin);
  fgets(buffer, 100, fin);

  // copy the needed number of children of the complete node
  for (int i = 0; i < nScenarios; ++i) {
    fgets(buffer, 100, fin);
    sscanf(buffer, "%d %d", &scen, &src);
    clusters[i] = src;

    // this node should be kept
    if (scen == src) {

      ttt = cNode->getChild(i);
      child = new Node(ttt->name());
      child->copy(ttt);
      rParent->addChild(child);
      keepnodes[i] = child;
    }
  }

  fclose(fin);

  // store the mapping of nodes
  for (int i = 0; i < nScenarios; ++i) {
    ttt = cNode->getChild(i);
    nMap[ttt] = keepnodes[clusters[i]];
  }

  // clean up
  delete[] clusters;

  return 0;
}

/**
 *  Generate a reduced tree by choosing a subset of scenarios.
 *
 *  @param nScenarios:
 *         Number of scenarios to appear in the reduced tree.
 *  @param clusteringFile:
 *         Name of the clustering file for scenario reduction.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::reduceTree(const int nScenarios, const char *clusteringFile) {

  printf(" --------------- reduceTree ----------------\n");

  // check that the number of scenarios is sensible
  if (nScenarios < 1) {
    printf("The number of scenarios must be positive.\n");
    return 1;
  }

  const Node *cNode = smps.getRootNode();

  // do not request more scenarios than there are in the complete tree
  int nWanted = nScenarios;
  if (nWanted > smps.getMaxScens())
    nWanted = smps.getMaxScens();

  // allocate the root node for the reduced tree
  rTree.setRootNode(new Node(100 + cNode->name()), cNode);
  Node *rNode = rTree.getRootNode();

  // copy the root node
  rNode->copy(cNode);
  nMap[cNode] = rNode;

  // print the number of scenarios in the reduced tree
  printf("Choosing %d scenarios.\n", nWanted);

  int rv = 0;

  // build up the reduced tree by selecting some scenarios
  if (clusteringFile)
    rv = reduceScenariosCluster(clusteringFile);
  else
    reduceScenarios(cNode, rNode, nWanted);

  if (rv)
    return rv;

  // recompute the probabilities in the reduced tree
  do {

    // find the corresponding node in the reduced tree
    rNode = nMap[cNode];

    // update its path probability by summing the path probabilities of
    // the complete nodes that map to it
    rNode->setProb(rNode->probNode() + cNode->probNode());

  } while ((cNode = cNode->next()));

  return 0;
}

/**
 *  Generate a reduced tree with aggregation.
 *
 *  @param nAggr:
 *         Number of stages to be aggregated (starting from the last).
 *  @return 1 If something goes wrong, 0 otherwise.
 */
int SmpsOops::aggregateStages(const int nAggr) {

  printf(" --------------- aggregateStages -----------\n");

  // check that the number of stages to aggregate is sensible
  const int last = smps.getPeriods() - nAggr + 1;
  if (nAggr < 2 || last < 2) {
    printf("No aggregation possible.\n");
    return 1;
  }

  // aggregate the last stages
  printf("Aggregating stages %d to %d.\n", last, smps.getPeriods());

  int nodeName = 100;
  const Node *cNode = smps.getRootNode();
  Node *rNode = NULL;

  // queue of nodes to be processed
  queue<const Node*> cNodes;
  cNodes.push(cNode);

  // copy the tree up to the aggregation point
  while (!cNodes.empty()) {

    // take the first element in the queue
    cNode = cNodes.front();
    cNodes.pop();
    assert(cNode != NULL);

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
  rTree.setRootNode(rNode, smps.getRootNode());

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
 *  Generate a subtree rooted at the given node.
 *
 *  @param cOrig:
 *         The node to be used as root node of the subtree.
 *  @param nodeName:
 *         A number used to differentiate the names of the subtree nodes.
 *  @return 1 If something goes wrong, 0 otherwise.
 */
int SmpsOops::createSubtree(const Node *cOrig, const int nodeName) {

  Node *rNode = NULL;
  const double rootProb = cOrig->probNode();

  // queue of nodes to be processed
  queue<const Node*> qNodes;
  qNodes.push(cOrig);

  // copy the nodes in the subtree
  while (!qNodes.empty()) {

    // take the first element in the queue
    const Node *cNode = qNodes.front();
    qNodes.pop();
    assert(cNode != NULL);

    // create a node in the reduced tree
    rNode = new Node(nodeName + cNode->name());
    rNode->copy(cNode);
    if (cNode != cOrig)
      nMap[cNode->parent()]->addChild(rNode);

    // scale the probability of the reduced node by the probability of
    // the root node of the subtree, so that the sum of the probabilities
    // of all nodes at the same stage is 1.0
    rNode->setProb(cNode->probNode() / rootProb);
    assert(rNode->probNode() <= 1.0);

#ifdef DEBUG_RTREE
    printf("Node %d: updating probability to %f\n",
           rNode->name(), rNode->probNode());
#endif

    // map the complete node to the reduced node
    nMap[cNode] = rNode;

    // put the children in the queue
    for (int i = 0; i < cNode->nChildren(); ++i)
      qNodes.push(cNode->getChild(i));
  }

  // find the root of the reduced tree
  while (rNode->parent())
    rNode = rNode->parent();

  // set the root of the reduced tree
  rTree.setRootNode(rNode, cOrig);

  return 0;
}

/**
 *  Set up the Oops algebras and vectors and build the primal-dual problem.
 *
 *  @param Pb:
 *         The SmpsReturn structure of the problem generated.
 *  @param wsIterate:
 *         The WSPoint to be used to warmstart the problem.
 *  @return Pointer to the problem to be solved by Oops.
 *
 *  @note
 *  InitAlgebrasNew assumes that a callback function is set up for all
 *  SparseMatrix leaves of A and Q.
 */
PDProblem SmpsOops::setupProblem(SmpsReturn &Pb, const WSPoint *wsIterate) {

  Algebra *A = Pb.AlgA;
  Algebra *Q = Pb.AlgQ;

  Algebra *AlgAug = InitAlgebrasNew(A, Q);

  Vector *vb = new Vector(A->Trow, "vb");
  Vector *vc = new Vector(A->Tcol, "vc");
  Vector *vu = new Vector(A->Tcol, "vu");
  Vector *vl = new Vector(A->Tcol, "vl");

  Vector *vx = new Vector(A->Tcol, "vx");
  Vector *vy = new Vector(A->Trow, "vy");
  Vector *vz = new Vector(A->Tcol, "vz");
  Vector *vs = NULL, *vw = NULL;
  if (smps.hasUpperBounds()) {
    vs = new Vector(A->Tcol, "vs");
    vw = new Vector(A->Tcol, "vw");
  }

  // use the warmstart point if available
  if (wsReady) {
    SmpsDenseToVector(wsIterate->x, vx, Pb, ORDER_COL);
    SmpsDenseToVector(wsIterate->y, vy, Pb, ORDER_ROW);
    SmpsDenseToVector(wsIterate->z, vz, Pb, ORDER_COL);
    if (smps.hasUpperBounds()) {
      SmpsDenseToVector(wsIterate->s, vs, Pb, ORDER_COL);
      SmpsDenseToVector(wsIterate->w, vw, Pb, ORDER_COL);
    }
  }

  vb->copyFromDense(Pb.b);
  vc->copyFromDense(Pb.c);
  vu->copyFromDense(Pb.u);
  vl->copyFromDense(Pb.l);

  // create the primal dual problem
  PDProblem pdProb(AlgAug, vb, vc, vu, vx, vy, vz);
  if (smps.hasUpperBounds()) {
    pdProb.s = vs;
    pdProb.w = vw;
  }
  pdProb.l = vl;

  return pdProb;
}

/**
 *  Set up a warmstart point from a reduced-tree solution.
 *
 *  This function takes the Vectors from the pdProb that contain the solution
 *  to the (reduced) problem just solved to setup the DenseVectors in wsPoint
 *  (for the complete problem).
 *
 *  This considers that the Vectors in pdProb contain the reshuffling done
 *  by generateSmps() to minimise the size of the Schur complement: this is
 *  removed by the calls to VectorToSmpsDense(), which produces DenseVectors
 *  of the dimension of the reduced problem. With the calls to setupWSPoint(),
 *  these are then expanded to the dimension of the complete tree.
 *
 *  @param pdProb:
 *         The reduced primal-dual problem.
 *  @param Ret:
 *         The SmpsReturn structure of the reduced problem.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::setupWarmStart(const PDProblem &pdProb, const SmpsReturn &Ret) {

  printf(" --------------- setupWarmStart ------------\n");

  // dimensions of the complete and the reduced deterministic equivalents
  const int nRows = smps.getTotRows(), rRows = rTree.getTotRows();
  const int nCols = smps.getTotCols(), rCols = rTree.getTotCols();

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
  WSPoint wsReduced(rRows, rCols, smps.hasUpperBounds());

  // recover the initial ordering of the solution vectors
  VectorToSmpsDense(pdProb.x, wsReduced.x, Ret, ORDER_COL);
  VectorToSmpsDense(pdProb.y, wsReduced.y, Ret, ORDER_ROW);
  VectorToSmpsDense(pdProb.z, wsReduced.z, Ret, ORDER_COL);
  if (smps.hasUpperBounds()) {
    VectorToSmpsDense(pdProb.s, wsReduced.s, Ret, ORDER_COL);
    VectorToSmpsDense(pdProb.w, wsReduced.w, Ret, ORDER_COL);
  }

  // allocate space for the vectors in the complete iterate
  if (!wsPoint) {
    wsPoint = new WSPoint(nRows, nCols, smps.hasUpperBounds());
  }

  // copy the solution from wsReduced into wsPoint
  int rv = setupWSPoint(wsPoint, &wsReduced, rTree.getOrigNode());
  if (rv)
    return rv;

  // erase the mapping
  nMap.clear();

  return 0;
}

/**
 *  Copy an available solution into a warmstart point.
 *
 *  This function traverses all nodes in the complete tree: for each of them,
 *  it finds the corresponding reduced tree node and copies its part of the
 *  solution, taking care of reconciling the probabilities between the two
 *  trees.
 *
 *  This function also works the other way round, that is it can initialise a
 *  reduced point from a complete one. In this case, as the mapping of nodes
 *  (nMap) always refers to nodes in the complete tree, we have to swap the
 *  roles of cNode and rNode.
 *
 *  @param wsNew:
 *         The WSPoint to be set up.
 *  @param wsRed:
 *         The WSPoint which provides the most current solution.
 *  @param rootNode:
 *         Node in the complete tree corresponding to the reduced tree root.
 *  @param swap:
 *         True if we copy the solution from a complete warmstart point.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsOops::setupWSPoint(const WSPoint *wsNew, WSPoint *wsRed,
                           const Node *rootNode, const bool swap) {

  const Node *cNode = rootNode, *rNode;
  int cIndex, rIndex, nElems;
  double crProb;

  double *xnew = wsNew->x->elts, *ynew = wsNew->y->elts, *znew = wsNew->z->elts;
  double *snew = wsNew->s ? wsNew->s->elts : NULL;
  double *wnew = wsNew->w ? wsNew->w->elts : NULL;
  double *xred = wsRed->x->elts, *yred = wsRed->y->elts, *zred = wsRed->z->elts;
  double *sred = wsRed->s ? wsRed->s->elts : NULL;
  double *wred = wsRed->w ? wsRed->w->elts : NULL;

  // ensure that we do not use the swap in the wrong case
  assert((wsNew == wsPoint) && (swap == false) ||
         (wsRed == wsPoint) && (swap == true));

  // go through the nodes in the complete tree
  do {

    rNode = nMap[cNode];
    if (!rNode)
      continue;

    // probability ratio between the nodes in the two trees used to
    // scale correctly the variables of the dual constraints (z, w, y)
    crProb = cNode->probNode() / rNode->probNode();

    cIndex = cNode->firstCol();
    rIndex = rNode->firstCol();
    nElems = cNode->nCols();
    if (swap) {
      cIndex = rNode->firstCol();
      rIndex = cNode->firstCol();
      crProb = rNode->probNode() / cNode->probNode();
    }

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
      xnew[cIndex + i] = xred[rIndex + i];
      znew[cIndex + i] = zred[rIndex + i] * crProb;
      if (smps.hasUpperBounds()) {
	snew[cIndex + i] = sred[rIndex + i];
	wnew[cIndex + i] = wred[rIndex + i] * crProb;
      }
    }

    cIndex = cNode->firstRow();
    rIndex = rNode->firstRow();
    nElems = cNode->nRows();
    if (swap) {
      cIndex = rNode->firstRow();
      rIndex = cNode->firstRow();
    }

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
      ynew[cIndex + i] = yred[rIndex + i] * crProb;
    }

  } while ((cNode = cNode->next()));

  return 0;
}

/**
 *  Retrieve the solution.
 *
 *  @param pdProb:
 *         Pointer to the problem solved by OOPS.
 *  @param Ret:
 *         The SmpsReturn structure of the problem.
 *
 *  @bug
 *  The value of the slacks is always zero, as the slacks would need to be
 *  computed here.
 */
int SmpsOops::getSolution(PDProblem &pdProb, SmpsReturn &Ret) {

  DenseVector *x, *y, *z, *r;
  Tree *Trow = Ret.AlgA->Trow;
  Tree *Tcol = Ret.AlgA->Tcol;
  const int nRows = Trow->end - Trow->begin;
  const int nCols = Tcol->end - Tcol->begin;

  // allocate space for the solution vectors
  x = NewDenseVector(nCols, "dx");
  y = NewDenseVector(nRows, "dy");
  z = NewDenseVector(nCols, "dz");

  // this is just a placeholder, we cannot yet compute the slacks
  r = NewDenseVector(nRows, "slacks");

  // recover the original ordering in the solution vectors
  VectorToSmpsDense(pdProb.x, x, Ret, ORDER_COL);
  VectorToSmpsDense(pdProb.y, y, Ret, ORDER_ROW);
  VectorToSmpsDense(pdProb.z, z, Ret, ORDER_COL);

  // print the solution
#ifdef WITH_MPI
  if (IS_ROOT_PAR)
#endif
    printSolution(Ret.rootNode, x->elts, y->elts, r->elts, z->elts);

  // clean up
  FreeDenseVector(x);
  FreeDenseVector(y);
  FreeDenseVector(z);
  FreeDenseVector(r);

  return 0;
}

/**
 *  Order the nodes of the tree according to the cutoff level.
 *
 *  The nodes of the event tree are reordered according to the cutoff level:
 *  - breadth-first until the given level (so that these nodes come first);
 *  - depth-first afterwards (to give diagonal structure in rest).
 *
 *  The cutoff level has to be positive and smaller than the number of
 *  stages in the problem. Also, it is capped to MAX_CUTOFF.
 *
 *  @param Tree:
 *         The tree to be reordered.
 *  @return A nonzero value if something goes wrong; 0 otherwise.
 */
int SmpsOops::orderNodes(SmpsTree &Tree) {

  const int nPeriods = smps.getPeriods();
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
  int nBlocks = 0;

  // queue of nodes to be followed
  queue<Node*> qNodes;

  // queue of nodes in the new order
  queue<Node*> qOrder;

  // there's nothing to order if the node doesn't have children, which is
  // the case when running the decomposition on a two-stage problem
  if (node->nChildren() == 0)
    goto common;

  // shift the cutoff by the period of the root node, which may not be zero
  // in the decomposition case
  cutoff += node->level();

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
    dfsNode(qOrder, node, nBlocks, cutoff);
  }

  // update the next links

  // start from the root
  node = Tree.getRootNode();

  // restore the original cutoff level
  cutoff -= node->level();

  do {

    qOrder.pop();
    node->setNext(qOrder.empty() ? NULL : qOrder.front());

  } while ((node = node->next()));

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

 common:

  // store the number of diagonal blocks
  Tree.setBlocks(nBlocks);

  // reset the period starts
  smps.setNodeStarts(Tree);

  // dimensions of the deterministic equivalent
  printf("The deterministic equivalent matrix is %dx%d, %d nonzeros.\n",
	 Tree.getTotRows(), Tree.getTotCols(), smps.countNonzeros(Tree));

  return 0;
}

/** Perform a recursive depth-first ordering of the node and its children */
void dfsNode(queue<Node*> &qOrder, Node *node,
             int &nBlocks, const int cutoff) {

  // put the node in the reordered queue
  qOrder.push(node);

  // count the number of diagonal blocks
  if (node->level() == cutoff)
    ++nBlocks;

  node->setBlock(nBlocks);

  // traverse the children in depth-first order
  for (int i = 0; i < node->nChildren(); i++) {
    dfsNode(qOrder, node->getChild(i), nBlocks, cutoff);
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
    FreeAlgebraAlg(AlgA);
  }

  if (AlgQ) {
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
WSPoint::WSPoint(const int nRows, const int nCols, const bool hasUpperBounds) :
  x(NewDenseVector(nCols, "x_ws")),
  y(NewDenseVector(nRows, "y_ws")),
  z(NewDenseVector(nCols, "z_ws")),
  s(NULL),
  w(NULL) {
  if (hasUpperBounds) {
    s = NewDenseVector(nCols, "s_ws");
    w = NewDenseVector(nCols, "w_ws");
  }
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
  _useDecomposition(0),
  _warmstartDecomposition(0),
  _cutoffLevel(1),
  _clusteringFileName(NULL) {
}

/** Parse the command line options */
int OptionsOops::parse() {

  // add the specialised options
  Options::addOption("-w", "use a warmstart strategy with scenario reduction",
		     &_useReduction, true);
  Options::addOption("-c", "specify a clustering file for scenario reduction",
                     &_clusteringFileName);
  Options::addOption("-a", "use a warmstart strategy with stage aggregation",
		     &_useAggregation, true);
  Options::addOption("-k", "use a warmstart strategy with tree decomposition",
                     &_useDecomposition);
  Options::addOption("-x", "warmstart the subproblems for decomposition",
                     &_warmstartDecomposition);
  Options::addOption("-l", "cutoff level (for multistage programs)",
		     &_cutoffLevel, true);

  // parse the common options
  int rv = Options::parse();

  return rv;
}

/**
 *  Create the first stage block from the core file.
 *
 *  This builds a sparsematrix containing the linking block between 1st and 2nd
 *  stage. It contains the 1st stage variables that appear in 2nd stage
 *  constraints.
 *
 *  @return The product Tx to be subtracted from the right-hand side in the
 *          subproblems.
 *
 *  @note
 *  This assumes that the linking blocks are identical for all second-stage
 *  nodes, so that the contribution is computed just once.
 */
double* SmpsOops::firstStageContribution() {

  const SparseData data = smps.getSparseData();

  const int nRows = smps.getBegPeriodRow(2) - smps.getBegPeriodRow(1);
  const int nCols = smps.getBegPeriodCol(1) - smps.getBegPeriodCol(0);

  // the number of nonzeros is an upper approximation here, since it contains
  // the nonzeros appearing above this block and in the objective
  const int nNonz = data.clpnts[nCols];
  const int firstRow = smps.getBegPeriodRow(1);

  if (!wsPoint) {
    printf("No wsPoint available to obtain the first stage correction!\n");
    return NULL;
  }

  SparseSimpleMatrix *sparse = NewSparseMatrix(nRows, nCols, nNonz,
                                               "1stageblock");
  sparse->nb_el = 0;

#ifdef DEBUG_DECOMPOSITION
  {
    const Node *node = smps.getRootNode()->getChild(0);
    printf("First-stage block is (%dx%d)\n", nRows, nCols);
    printf("root node is node %d: rows %d, cols %d\n", node->name(),
           node->nRows(), node->nCols());
  }
#endif

  // copy the elements in the block
  for (int col = 0; col < nCols; ++col) {

    sparse->col_beg[col] = sparse->nb_el;

    // copy the nonzeros in the current column
    for (int el = data.clpnts[col]; el < data.clpnts[col + 1]; ++el) {

      // this element belongs to a first stage constraint
      if (data.rwnmbs[el] < firstRow) {
        /*
        printf("Skipping element %d: coeff %f (row %d)\n",
               el, data.acoeff[el], data.rwnmbs[el]);
        */
        continue;
      }

      assert(sparse->nb_el < sparse->max_nb_el);
      sparse->element[sparse->nb_el] = data.acoeff[el];
      sparse->row_nbs[sparse->nb_el] = data.rwnmbs[el] - firstRow;

#ifdef DEBUG_DECOMPOSITION
      printf("Col %d:: el: %d (max: %d) row: %d  coeff: %f\n", col,
             sparse->nb_el, sparse->max_nb_el, sparse->row_nbs[sparse->nb_el],
             sparse->element[sparse->nb_el]);
#endif

      assert(sparse->row_nbs[sparse->nb_el] >= 0);
      assert(sparse->row_nbs[sparse->nb_el] < sparse->nb_row);

      sparse->nb_el++;
      sparse->col_len[col]++;
    }
  }

  sparse->col_beg[nCols] = sparse->nb_el;

  Algebra *sp = NewSparseSimpleAlgebra(sparse);
  Tree tRow(0, nRows, 0);
  Tree tCol(0, nCols, 0);
  tRow.setLeavesLocal(NULL);
  tCol.setLeavesLocal(NULL);
  tRow.setIndex();
  tCol.setIndex();

  // get the first stage vector
  Vector *v = new Vector(&tCol, "v", wsPoint->x->elts); // XXX does the ordering here matter?
  Vector *vsol = new Vector(&tRow, "sol");

  // vsol = sp * v
  sp->MatrixTimesVect(sp, v, vsol, 0, 1.0);

#ifdef DEBUG_DECOMPOSITION
  PrintVector(v);
  sp->Print(stdout, sp->Matrix, "%f");
  PrintVector(vsol);
#endif

  DenseVector *ddense = GetDenseVectorFromVector(vsol);
  double *array = new double[nRows];
  memcpy(array, &ddense->elts[0], nRows * sizeof(double));

  // clean up
  FreeAlgebraAlg(sp);
  delete v;
  delete vsol;

  return array;
}
