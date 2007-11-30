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

#include "oops/oops.h"
#include "Smps.h"
#include "options.h"


// forward declaration
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

  /** Solve the problem instance */
  int solve(const opt_st *options, HopdmOptions *hopdm_options);

  int getSolution(primal_dual_pb *Prob, SmpsReturn *Ret);

 private:

  /** The smps instance to solve */
  Smps smps;

  /** The cutoff level */
  int level;

  /** Number of diagonal blocks in the deterministic equivalent */
  int nBlocks;

  int *block;

  int *order;

  int *revorder;

  /** Generate the deterministic equivalent for the smps instance */
  SmpsReturn* generateSmps(void);

  /** Return information on the nodes necessary for printing the solution */
  const NodeInfo* getNodeInfo(void) const;

  /** Order the nodes according to the cutoff level */
  void orderNodes(void);

  /** Apply the scenario changes */
  int applyScenarios(SmpsReturn *Ret,
		     Algebra **DiagEntries, Algebra **RightColEntries,
		     int *f_rw_blk, int *f_cl_blk);

  void reorderObjective(DenseVector *obj, DenseVector *upb,
			int *is_col_diag, const int rnkn, char **colnames);

  void addChildrenToList(const int node, int *next);

  void setNodeChildrenRnkc(Algebra **RC, Algebra **DG,
			   int *p_pd_rw, int *f_rw_blk, int *is_col_diag,
			   const SparseData &data,
			   const int node, const int colBlk,
			   const int rnkCol, const int coreCol);

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

  /** Rows of deterministic equivalent */
  int ttm;

  /** Columns of deterministic equivalent */
  int ttn;

  /** Columns in Rnk part of RankCor */
  int nb_col_rnk;

  /** Columns in D0 part of RnkCor */
  int nb_col_diag;

  /** Rows in first periods (making up [D0Rnk]) */
  int nb_row_rnk;

  /** Columns in first 'level' periods of core matrix that can be shifted
      to the first diagonal block */
  int *is_col_diag;

  /** Period array from SmpsTree structure */
  int *tree_period;

  /** Order array from SmpsTree structure */
  int *tree_order;

  /** f_cl_pd array from SmpsTime structure */
  int *time_f_cl_pd;

};


/** The orientation of a vector needed for reordering */
enum {
  ORDER_ROW,
  ORDER_COL
};

/** Free the space allocated for the SmpsReturn structure */
void freeSmpsReturn(SmpsReturn *ret);

void SmpsVectorToDense(Vector *x, DenseVector *dx,
		       SmpsReturn *Ret, const int rowcol);
void SmpsDenseToVector(DenseVector *dx, Vector *x,
		       SmpsReturn *Ret, const int rowcol);

#endif /* _SMPS_OOPS_H_ */
