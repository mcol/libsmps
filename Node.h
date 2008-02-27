/**
 *  @file Node.h
 *
 *  Definition of the Node object.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
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
  Node(const int name);

  /** Destructor */
  ~Node();

  /** Add a node to the list of children */
  int addChildNode(Node *child);

  /** Retrieve the child node at the given index */
  Node *getChild(const int childNumber) const;

  // Setters

  /** Set the scenario index */
  void setScen(const int scen) {
    _scenario = scen;
  }

  /** Set the path probability for this node */
  void setProb(const double prob) {
    _probNode = prob;
  }

  /** Set the indices to the matrix data */
  void setMatrixPointers(const int row, const int col,
			 const int nRows, const int nCols) {
    _firstRow = row;
    _firstCol = col;
    _nRows = nRows;
    _nCols = nCols;
  }

  // Accessor functions

  /** Return the name of the node */
  int name(void) const {
    return _name;
  }

  /** Return the number of children of the node */
  int nChildren(void) const {
    return _children.size();
  }

  /** Return the level of the node in the tree */
  int level(void) const {
    return _level;
  }

  /** Return the scenario index of the node */
  int scenario(void) const {
    return _scenario;
  }

  /** Return the probability associated to the node */
  double probNode(void) const {
    return _probNode;
  }

  /** Print some of the information of the node */
  void printNode(void) const;

 private:

  /** Pointer to the parent node */
  Node *_parent;

  /** Set of pointers to the children nodes */
  vector<Node*> _children;

  /** Index name of this node */
  int _name;

  /** Level of the node within a tree */
  int _level;

  /** Scenario this node belongs to */
  int _scenario;

  /** Path probability of this node */
  double _probNode;

  /** Index of the first row of node (in the deterministic equivalent) */
  int _firstRow;

  /** Index of the first column of node (in the deterministic equivalent) */
  int _firstCol;

  /** Number of rows for this node */
  int _nRows;

  /** Number of columns for this node */
  int _nCols;

};

#endif /* _NODE_H_ */
