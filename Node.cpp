/*
 *  Node.cpp
 *
 *  The Node object.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <vector>
#include <stdio.h>
#include <assert.h>
#include "Node.h"
using namespace std;


/** Constructor */
Node::Node(const int nodeName) :
  _parent(NULL),
  _next(NULL),
  _name(nodeName),
  _level(0),
  _nLevels(1),
  _scenario(0),
  _block(0),
  _probNode(0.0),
  _firstRow(0),
  _firstCol(0),
  _nRows(0),
  _nCols(0) {
}

/** Destructor */
Node::~Node() {

  // free the space occupied by the children nodes
  for (int i = 0; i < (int) _children.size(); ++i) {
    delete _children[i];
  }
}

/** Copy the information from a node */
int Node::copy(const Node *fromNode) {

  _level    = fromNode->_level;
  _scenario = fromNode->_scenario;

  return 0;
}

/**
 *  Check whether the node is part of the given scenario.
 *
 *  Each node is associated to a specific scenario, which is stored in
 *  the member variable _scenario. However, a non-leaf node belongs to
 *  all the scenarios that its children belong to, as it appears in all
 *  those scenarios.
 *
 *  @param scen:
 *         The index of the scenario we are interested in.
 *  @return True is the node belongs to the given scenario; False otherwise.
 */
bool Node::belongsToScenario(const int scen) const {

  // check the scenario of this node
  if (_scenario == scen)
    return true;

  // check the scenarios of the children
  for (int i = 0; i < (int) _children.size(); ++i)
    if (_children[i]->_scenario == scen)
      return true;

  return false;
}

/** Add a node to the list of children */
int Node::addChild(Node *childNode) {

  // ignore the request to add a child to a NULL node
  if (!this)
    return 1;

  // don't add a node as a child of itself
  assert(childNode != this);

  // link the child to the parent
  childNode->_level = _level + 1;
  childNode->_parent = this;

  // store the new child
  _children.push_back(childNode);

#ifdef DEBUG_NODE
  printf("Adding child %d to node %d\n", childNode->name(), this->name());
#endif

  return 0;
}

/** Remove a node from the list of children */
int Node::removeChild(Node *childNode) {

  for (int i = 0; i < nChildren(); ++i) {

    // compare the addresses
    if (_children[i] == childNode) {
      _children.erase(_children.begin() + i);
      return 0;
    }
  }

  // could not find the child
  return 1;
}

/** Retrieve the child node at the given index */
Node* Node::getChild(const int childNumber) const {

  // check that the index provided is sensible
  if (childNumber > (int) _children.size() - 1 ||
      childNumber < 0)
    return NULL;

  return _children[childNumber];
}

/** Print some of the information of the node */
void Node::print() const {

  // node parent scen n_chd per prob | rows cols
  printf(" %5i %5i %5i  %4i  %4i   %.4f  | %4i  %4i  (%4i)\n",
	 _name,
	 _parent ? _parent->name() : 0,
	 _scenario + 1,
	 nChildren(),
	 _level + 1,
	 _probNode,
	 _nRows,
	 _nCols,
	 next() ? next()->name() : -1);
}
