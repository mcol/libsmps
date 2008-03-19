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
#include <assert.h>
#include "Node.h"
using namespace std;


/** Constructor */
Node::Node(const int nodeName) :
  _parent(NULL),
  _next(NULL),
  _name(nodeName),
  _level(0),
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

/** Add a node to the list of children */
int Node::addChild(Node *childNode) {

  // don't add a node as a child of itself
  assert(childNode != this);

  // link the child to the parent
  childNode->_level = _level + 1;
  childNode->_parent = this;

  // store the new child
  _children.push_back(childNode);

  return 0;
}

/** Retrieve the child node at the given index */
Node* Node::getChild(const int childNumber) const {

  // check that the index provided is sensible
  if (childNumber > (int) _children.size() ||
      childNumber < 0)
    return NULL;

  return _children[childNumber];
}

/** Print some of the information of the node */
void Node::print() const {

  // node parent scen n_chd per prob | rows cols
  printf("  %4i  %4i  %4i  %4i  %4i   %.4f  | %4i  %4i  (%4i)\n",
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
