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

#include <queue>
#include "oops/oops.h"
#include "Smps.h"
#include "interface.h"


// forward declarations
class OptionsOops;
struct SmpsReturn;


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
  int solve(const OptionsOops &opt, HopdmOptions *hopdm_options);

  /** Solve a reduced problem */
  int solveReduced(const OptionsOops &opt, HopdmOptions *hopdm_options);

  int getSolution(PDProblem *Prob, SmpsReturn *Ret);

 private:

  /** The smps instance to solve */
  Smps smps;

  /** Root node of a reduced event tree */
  Node *rootReduced;

  /** The cutoff level */
  int level;

  /** Number of diagonal blocks in the deterministic equivalent */
  int nBlocks;

  /** Return the root node of the reduced event tree */
  Node* getRootNodeReduced(void) {
    return rootReduced;
  }

  /** Generate the deterministic equivalent for the smps instance */
  SmpsReturn* generateSmps(const Node *rootNode);

  /** Set up the Oops algebras and vectors and build the primal-dual problem */
  PDProblem* setupProblem(SmpsReturn *Pb);

  /** Order the nodes according to the cutoff level */
  int orderNodes(Node *node);

  /** Perform a recursive depth-first ordering of the node and its children */
  void dfsNode(queue<Node*> &qOrder, Node *node);

  /** Apply the scenario changes */
  int applyScenarios(SmpsReturn *Ret,
		     Algebra **DiagEntries, Algebra **RightColEntries,
		     int *f_rw_blk, int *f_cl_blk);

  /** Reorder objective and bounds */
  void reorderObjective(SmpsReturn *Ret, const int rnkn);

  void setNodeChildrenRnkc(Algebra **RC, Algebra **DG,
			   int *p_pd_rw, int *f_rw_blk, int *is_col_diag,
			   const SparseData &data,
			   const Node *node, const int colBlk,
			   const int rnkCol, const int coreCol);

  /** Copy a Vector into a breadth-first ordered DenseVector */
  void VectorToSmpsDense(Vector *v, DenseVector *dv,
			 const SmpsReturn *Ret, const int ordering);

  /** Copy a breadth-first ordered DenseVector into a Vector */
  void SmpsDenseToVector(DenseVector *dv, Vector *v,
			 const SmpsReturn *Ret, const int ordering);

  void backOrderRowVector(double *x, const SmpsReturn *Ret);
  void backOrderColVector(double *x, const SmpsReturn *Ret);
  void forwOrderRowVector(double *x, const SmpsReturn *Ret);
  void forwOrderColVector(double *x, const SmpsReturn *Ret);

  /** Free the space allocated for the SmpsReturn structure */
  void freeSmpsReturn(SmpsReturn *ret);

};


/**
 *  The returned data from generateSmps.
 *
 *  This structure is used to pass the Algebras and Vectors from the
 *  SMPS generator to OOPS.
 */
struct SmpsReturn {

  /** Algebra for the linear part */
  Algebra *AlgA;

  /** Algebra for the quadratic part */
  Algebra *AlgQ;

  /** Right-hand side vector */
  DenseVector *b;

  /** Objective vector */
  DenseVector *c;

   /** Upper bounds vector */
  DenseVector *u;

  /** Names of the rows in the deterministic equivalent */
  char **rownames;

  /** Names of the columns in the deterministic equivalent */
  char **colnames;

  /* The following entries are for the ordering/re-ordering routines
     (these values have to be remembered from the generating phase) */

  /** Columns in Rnk part of RankCor */
  int nb_col_rnk;

  /** Columns in D0 part of RnkCor */
  int nb_col_diag;

  /** Rows in first periods (making up [D0Rnk]) */
  int nb_row_rnk;

  /** Columns in first 'level' periods of core matrix that can be shifted
      to the first diagonal block */
  int *is_col_diag;

  /** Root node of the tree for which the problem has been generated */
  const Node *rootNode;

};


/** Command line options for the Oops interface */
class OptionsOops : public Options {

  /** Whether a warmstart strategy should be used */
  int _useWarmstart;

  /** The value of the cutoff level (for multistage problems) */
  int _cutoffLevel;

 public:

  /** Constructor */
  OptionsOops(const int argc = 0, const char *argv[] = NULL);

  /** Parse the command line options */
  int parse(void);

  /** Retrieve the value of the useWarmstart option */
  int useWarmstart(void) const {
    return _useWarmstart;
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
