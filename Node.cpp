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
#include "Node.h"
using namespace std;


/** Constructor */
Node::Node() :
  level(0),
  nChildren(0),
  parent(NULL) {
}

/** Destructor */
Node::~Node() {

  // free the space occupied by the children nodes
  for (vector<Node *>::iterator childNode = children.begin();
       childNode != children.end(); ++childNode) {
    delete *childNode;
  }
}

/** Add a node to the list of children */
int Node::addChild() {

  Node *child = new Node;
  child->level = level + 1;
  child->parent = this;
  children.push_back(child);

  return 0;
}
