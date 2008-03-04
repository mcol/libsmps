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

  // Setters

  /** Set the link to the next node in order */
  void setNext(Node *nextNode) {
    _next = nextNode;
  }

  /** Set the scenario index */
  void setScen(const int scen) {
    _scenario = scen;
  }

  /** Set the block number */
  void setBlock(const int blk) {
    _block = blk;
  }

  /** Set the path probability for this node */
  void setProb(const double prob) {
    _probNode = prob;
  }

  /** Set the indices to the matrix data */
  void setMatrixPointers(const int row, const int col,
			 const int rows, const int cols) {
    _firstRow = row;
    _firstCol = col;
    _nRows = rows;
    _nCols = cols;
  }

  // Accessor functions

  /** Return the parent node */
  Node* parent() const {
    return _parent;
  }

  /** Fast access to the next node in order */
  Node* next() const {
    return _next;
  }

  /** Retrieve the child node at the given index */
  Node* getChild(const int childNumber) const;

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

  /** Return the block number of the node */
  int block(void) const {
    return _block;
  }

  /** Return the probability associated to the node */
  double probNode(void) const {
    return _probNode;
  }

  /** Return the first row associated with this node */
  int firstRow(void) const {
    return _firstRow;
  }

  /** Return the first column associated with this node */
  int firstCol(void) const {
    return _firstCol;
  }

  /** Return the number of rows associated with this node */
  int nRows(void) const {
    return _nRows;
  }

  /** Return the number of columns associated with this node */
  int nCols(void) const {
    return _nCols;
  }

  /** Print some of the information of the node */
  void print(void) const;

 private:

  /** Pointer to the parent node */
  Node *_parent;

  /** Pointer to the next node in order */
  Node *_next;

  /** Set of pointers to the children nodes */
  vector<Node*> _children;

  /** Index name of this node */
  int _name;

  /** Level of the node within a tree */
  int _level;

  /** Scenario this node belongs to */
  int _scenario;

  /** Diagonal block this node belongs to */
  int _block;

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
