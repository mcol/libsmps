/*
 *  Node.h
 *
 *  Definition of the Node object.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#ifndef _NODE_H_
#define _NODE_H_

#include <vector>
using namespace std;


/** The Node class */
class Node {

 public:

  /** Constructor */
  Node();

  /** Destructor */
  ~Node();

  /** Add a node to the list of children */
  int addChild();

  /** Set the parent node */
  void setParent(Node *parentNode) {
    parent = parentNode;
  }

 private:

  /** The level of the node within a tree */
  int level;

  /** The number of children of the node */
  int nChildren;

  /** Pointer to the parent node */
  Node *parent;

  /** The set of pointers to the children nodes */
  vector<Node*> children;

};

#endif /* _NODE_H_ */
