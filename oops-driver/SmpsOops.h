/**
 *  @file SmpsOops.h
 *
 *  Definitions for the SMPS interface to Oops.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
 *
 */

#ifndef _SMPS_OOPS_H_
#define _SMPS_OOPS_H_

#include <map>
#include <queue>
#include "oops/oops.h"
#include "Smps.h"
#include "interface.h"


// forward declarations
class OptionsOops;
class WSPoint;
class SmpsReturn;

// the maximum cutoff level
// note that the current implementation doesn't exploit a cutoff larger than 1
static const int MAX_CUTOFF = 2;


/** The SmpsOops class */
class SmpsOops {

 public:

  /** Constructor */
  SmpsOops(string smpsFile, const int level = 1);

  /** Destructor */
  ~SmpsOops(void);

  /** Read the smps files */
  int read(void);

  /** Solve the complete problem */
  int solve(const OptionsOops &opt, HopdmOptions &hopdmOpts);

  /** Solve a reduced problem */
  int solveReduced(const OptionsOops &opt, HopdmOptions &hopdmOpts);

  /** Generate a reduced tree by choosing a subset of scenarios */
  int reduceTree(const int nScenarios);

  /** Generate a reduced tree with aggregation */
  int aggregateStages(const int nAggr);

  /** Retrieve the solution information */
  int getSolution(PDProblem *Prob, SmpsReturn &Ret);

 private:

  /** The smps instance to solve */
  Smps smps;

  /** The reduced event tree */
  SmpsTree rTree;

  /** A warmstart point for the complete problem */
  WSPoint *wsPoint;

  /** The mapping between the nodes in the complete and reduced trees */
  map<const Node*, Node*> nMap;

  /** The cutoff level */
  int cutoff;

  /** Number of diagonal blocks in the deterministic equivalent */
  int nBlocks;

  /** Generate the deterministic equivalent for the smps instance */
  int generateSmps(const SmpsTree &tree, SmpsReturn &ret);

  /** Set up the Oops algebras and vectors and build the primal-dual problem */
  PDProblem* setupProblem(SmpsReturn &Pb);

  /** Create a reduced tree in a recursive manner */
  void reduceScenarios(const Node *cNode, Node *rParent, const int nWanted);

  /** Set up the warmstart point from a reduced-tree solution */
  int setupWarmStart(const PDProblem *pdProb, const SmpsReturn &Pb);

  /** Recompute the probabilities in the reduced tree */
  void adjustProbabilities(void);

  /** Order the nodes according to the cutoff level */
  int orderNodes(SmpsTree &Tree);

  /** Perform a recursive depth-first ordering of the node and its children */
  void dfsNode(queue<Node*> &qOrder, Node *node);

  /** Apply the scenario changes */
  int applyScenarios(const SmpsTree &tree, SmpsReturn *Ret,
                     Algebra **Diagon, Algebra **Border,
		     int *f_rw_blk, int *f_cl_blk);

  /** Reorder objective and bounds */
  void reorderObjective(const SmpsTree &tree, SmpsReturn *Ret, const int rnkn);

  /** Copy a Vector into a breadth-first ordered DenseVector */
  void VectorToSmpsDense(Vector *v, DenseVector *dv,
                         const SmpsReturn &Ret, const int ordering);

  /** Copy a breadth-first ordered DenseVector into a Vector */
  void SmpsDenseToVector(DenseVector *dv, Vector *v,
                         const SmpsReturn &Ret, const int ordering);

};


/**
 *  The returned data from generateSmps.
 *
 *  This structure is used to pass the Algebras and Vectors from the
 *  SMPS generator to OOPS.
 */
class SmpsReturn {

 public:

  /** Constructor */
  SmpsReturn();

  /** Destructor */
  ~SmpsReturn();

  /** Algebra for the linear part */
  Algebra *AlgA;

  /** Algebra for the quadratic part */
  Algebra *AlgQ;

  /** Right-hand side vector */
  DenseVector *b;

  /** Objective vector */
  DenseVector *c;

  /** Lower bounds vector */
  DenseVector *l;

  /** Upper bounds vector */
  DenseVector *u;

  /** Names of the rows in the deterministic equivalent */
  char **rownames;

  /** Names of the columns in the deterministic equivalent */
  char **colnames;

  /* The following entries are for the ordering/re-ordering routines
     (these values have to be remembered from the generating phase) */

  /** Root node of the tree for which the problem has been generated */
  const Node *rootNode;

  /** Columns in first 'level' periods of core matrix that can be shifted
      to the first diagonal block */
  int *is_col_diag;

  /** Rows in first periods (making up [D0Rnk]) */
  int nRowsRnkc;

  /** Columns in Rnk part of RankCor */
  int nColsRnkc;

  /** Columns in D0 part of RnkCor */
  int nColsDiag;

};


/** A warmstart point for the complete problem */
class WSPoint {

 public:

  /** Constructor */
  WSPoint(DenseVector *vx, DenseVector *vy, DenseVector *vz,
	  DenseVector *vs = NULL, DenseVector *vw = NULL);

  /** Destructor */
  ~WSPoint();

  /** Primal vector */
  DenseVector *x;

  /** Dual vector */
  DenseVector *y;

  /** Dual slack vector */
  DenseVector *z;

  /** Upper bound slack vector */
  DenseVector *s;

  /** Dual upper bound vector */
  DenseVector *w;

 };


/** Command line options for the Oops interface */
class OptionsOops : public Options {

  /** Whether scenario reduction should be used */
  int _useReduction;

  /** Whether stage aggregation should be performed */
  int _useAggregation;

  /** The value of the cutoff level (for multistage problems) */
  int _cutoffLevel;

 public:

  /** Constructor */
  OptionsOops(const int argc = 0, const char *argv[] = NULL);

  /** Parse the command line options */
  int parse(void);

  /** Determine whether a warmstart strategy has to be employed */
  int useWarmstart(void) const {
    return _useReduction || _useAggregation;
  }

  /** Retrieve the value of the useReduction option */
  int useReduction(void) const {
    return _useReduction;
  }

  /** Retrieve the value of the useAggregation option */
  int useAggregation(void) const {
    return _useAggregation;
  }

  /** Retrieve the value of the cutoffValue option */
  int cutoffLevel(void) const {
    return _cutoffLevel;
  }

};


/** The orientation of a vector needed for reordering */
enum {
  ORDER_ROW,
  ORDER_COL
};

#endif /* _SMPS_OOPS_H_ */
