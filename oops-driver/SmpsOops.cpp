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
  orderNodes();

  return rv;
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
void SmpsOops::orderNodes() {

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
