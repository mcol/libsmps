/*
 *  generateSmps.cpp
 *
 *  Generate the Oops structures for the SMPS problem.
 *
 *  Andreas Grothey
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <string.h>
#include <assert.h>
#include "SmpsOops.h"


static void
setupRhs(const Smps &smps, SmpsReturn *Ret);

static void
setupObjective(const Smps &smps, SmpsReturn *Ret);

static void
backOrderColVector(const Smps &smps, const SmpsReturn &Ret, double *x);

static void
forwOrderColVector(const Smps &smps, const SmpsReturn &Ret, double *x);

static void
backOrderRowVector(const SmpsReturn &Ret, double *x);

static void
forwOrderRowVector(const SmpsReturn &Ret, double *x);

static int
copyLinkingBlocks(Smps &smps, Algebra **Array, const SparseData &data,
                  const Node *node, int *f_rw_blk,
                  int &idxCol, int &cIndex, int &nnzCol);

/**
 *  Generate the OOPS structures for the SMPS problem.
 *
 *  Uses the information from the SMPS files to create the deterministic
 *  equivalent problem.
 *
 *  All nodes in periods < cutoff are treated in a RankCorrector, nodes
 *  in period == cutoff are seeds for diagonal blocks.
 *
 *  Some extra analysis is done to:
 *  - remove columns from RankCor that do not affect any periods >= cutoff,
 *    these can be (and are) treated in an additional diagonal block (D-0)
 *    FIXME: should also separate rows from the first block of RankCor that
 *           do not affect columns in D-0
 *
 *  The subroutine procedes as follows:
 *
 *  (1a) Analyse which columns of the first block (periods < cutoff) can be
 *       treated in additional Diagonal D-0 and which remain in RnkC
 *       Count dimensions of D-0 and RnkC
 *	 Count Nonzeros in resulting RnkC and D-0 parts of big matrix
 *  (1b) Count Total nonzeros in big matrix
 *       Count Dimension & nonzeros of all RnkC and Diagonal parts
 *	 Set start row/col of Diagonal parts in big Matrix.
 *
 *  (2a) Allocate all Vectors and SparseMatrices
 *  (2b) Fill in correct bits of Big Matrix column by column
 *  (2c) Create big RHS, Obj, UB vectors (in order [blk-0 blk-1, ..., blk-n])
 *
 *  ( 3) Apply scenarios corrections to RHS, Obj and big Matrix
 *
 *  (4b) reorder Obj, UB vectors due to splitting blk-0 into
 *       D-0 and RnkC and using new order
 *                [ D-0 blk-1 ... blk-n RnkC]
 *
 *  (5a) Create Big Matrix out of SparseMatrices as RankCorAlgebra
 *  (4c) Correct RHS and objvalue.
 *	 Do this on whole Matrix (therefore needs to be done after 5a)
 */
int SmpsOops::generateSmps(const SmpsTree &tree, SmpsReturn &Ret) {

  // extract the root node
  const Node *rootNode = tree.getRootNode();

  // ensure that the root node exists
  if (!rootNode)
    return 1;

  int i, j;
  const Node *node;

  // store the root node
  Ret.rootNode = rootNode;

  // indicates columns in first periods that should be in RankCor part
  // but are not used outside these periods, so they are taken out of the
  // RankCor and form an additional diagonal block D-0
  int *is_col_diag;

  // nonzeros and number of columns in the new diagonal block D-0
  int nzdg0 = 0, cldg0 = 0;

  // period of the given node
  int perNode;

  // block of the given node
  int blkNode;

  // number of periods
  const int nPeriods = smps.getPeriods();

  // dimensions and nonzeros of rank corrector blocks
  int *rnkc_m_blk  = new int[nBlocks + 1];
  int *rnkc_n_blk  = new int[nBlocks + 1];
  int *rnkc_nz_blk = new int[nBlocks + 1];

  // dimensions and nonzeros of diagonal blocks
  int *diag_m_blk  = new int[nBlocks + 1];
  int *diag_n_blk  = new int[nBlocks + 1];
  int *diag_nz_blk = new int[nBlocks + 1];

  // initialise to zero
  memset(rnkc_nz_blk, 0, (nBlocks + 1) * sizeof(int));
  memset(diag_m_blk,  0, (nBlocks + 1) * sizeof(int));
  memset(diag_n_blk,  0, (nBlocks + 1) * sizeof(int));
  memset(diag_nz_blk, 0, (nBlocks + 1) * sizeof(int));

  // dimensions of the deterministic equivalent
  const int ttm = tree.getTotRows();
  const int ttn = tree.getTotCols();

  // first row/col in the diagonal part of the deterministic equivalent
  const int f_rwdiag = smps.getBegPeriodRow(cutoff);
  const int f_cldiag = smps.getBegPeriodCol(cutoff);

  // index of the objective row
  const int objRow = smps.getObjRowIndex();

  const SparseData data = smps.getSparseData();

  printf(" --------------- generateSmps --------------\n");

  /*
     Columns within periods <= cutoff are designated to go into the border
     part of the matrix. However some of them have no entries in rows
     corresponding to other time periods, and they can be treated differently.
     This code identifies these columns (is_col_diag[i]) and counts them.
     These columns can then be considered part of a new diagonal block.

     Since periods other than the first one may be in the border and non-first
     period columns of the core matrix will be repeated according to the
     scenario tree, we will need to count how many border columns go to the
     diagonal from every period.

     Scan through the columns of RankCorrector block to see which can go in
     the diagonal.
     For all row periods <= cutoff, count the number of nonzeros in columns
     that can go in the diagonal part (used to obtain the number of nonzeros
     in diagonal entry).
     For all col periods <= cutoff, count the number of columns that
     can go in the diagonal part (used to obtain col of diag).

     nzpddg[0-(cutoff-1)]  number of nonzeros in above columns
                           (this time sorted by their row-period in CORE)
                          FIXME: is this correct? Do we not need to
                                 distiguish by both row/col period
  */
  {
    is_col_diag = Ret.is_col_diag = new int[f_cldiag];

    // nonzeros in this column for each row-period
    int *nzpdd0 = new int[nPeriods];

    // total number of columns that have no entries in rows associated with
    // diagonal blocks and so can go into a new diagonal block
    int clpddg[MAX_CUTOFF];

    // total number of entries in the above columns for each row-period
    int nzpddg[MAX_CUTOFF];

    // whether the column has entries in rows associated with diagonal blocks
    bool found;

    int row;

    // initialise the vectors
    for (j = 0; j < cutoff; ++j)
      nzpddg[j] = clpddg[j] = 0;

    // for all columns in border block
    for (i = 0; i < f_cldiag; ++i) {

      found = false;
      is_col_diag[i] = 0;

      // count the number of nonzero entries in the column and their
      // distribution into row period blocks

      for (j = 0; j < cutoff; ++j)
	nzpdd0[j] = 0;

      // for all elements in this column
      for (j = data.clpnts[i]; j < data.clpnts[i + 1]; ++j) {

	row = data.rwnmbs[j];

	// not objective
	if (row != objRow) {

	  // count nonzeros this column will create in final matrix
	  nzpdd0[smps.getRowPeriod(row)]++;
	}

	// this entry is associated with a diagonal block
	if (row >= f_rwdiag)
	  found = true;
      }

      // no linking elements have been found
      if (!found) {

	// the column can go in the diagonal part
	is_col_diag[i] = 1;
	clpddg[smps.getColPeriod(i)]++;

	for (j = 0; j < cutoff; ++j)
	  nzpddg[j] += nzpdd0[j];
      }
    }

    // scan through all nodes up to the cutoff level and count the total
    // number of columns that can be moved into a diagonal block and the
    // corresponding number of nonzero elements

    node = rootNode;

    do {

      perNode = node->level();
      if (perNode >= cutoff)
	break;

#ifdef DEBUG_GENERATE_SMPS
      printf("Node %d period %d,  nz in diag %d, cols in diag %d\n",
	     node->name(), perNode, nzpddg[perNode], clpddg[perNode]);
#endif
      nzdg0 += nzpddg[perNode];
      cldg0 += clpddg[perNode];

    } while (node = node->next());

#ifdef DEBUG
    printf("Moved %d border columns (%d nonzeros) into a diagonal block.\n",
	   cldg0, nzdg0);
#endif

    // clean up
    delete[] nzpdd0;
  }

  /* Predict sizes of components of big Matrix:
     The matrix generated by DblBordDiagAlgebra looks like this

     [ D0            C0   ]
     [    D1         C1   ]
     [                    ]
     [          Dn-1 Cn-1 ]

     However due to the interpretation of parts:
      D1-Dn-1 : proper diagonal parts
      C1-Cn-1 : proper RankCorrector
      D0      : first periods part moved to diagonal
      C0      : first periods part treated as RankCor

     It is easier to think of it initially as a reordering

     [ D1               C1  ]
     [    D2            C2  ]
     [                   :  ]
     [          Dn-1    Cn-1]
     [               D0 C0  ]

     where [D0 C0] is the first periods matrix reordered to take into
     account columns that can be treated separately.

     For setting up the matrix the ordering does not make any difference at
     all, however when setting up RHS, Obj, UB vectors the ordering
     must be taken into account.

     Column vectors are set up in the order [ D0 D1 D2 ... Dn-1 Cn-1 ]
     and row vectors in the order [ C0 C1 C2 ... Cn-1 ]'
     (as in the first of the above matrices)

     Therefore the rows (and row vector b) are set up in the correct order,
     whereas columns (and column vectors c, obj, u) need to be reordered
     (and ordered back after the solution).

     B1 - Bn are diagonal blocks, A1 - An corresponding Border Blocks

     [ C00             ]
     [ C01 C11         ]
     [ C02 C12 C22     ]
     [ C03 C13 C23 C33 ]

                           root
                       /          \
                      /            \
                     /              \
                   1                  2
                /    \              /    \
               /      \            /      \
              3        4          5        6
            /   \    /   \      /   \    /   \
           7     8  9     10  11    12  13   14

   And cutoff level = 2 then the resulting matrix looks like this

       [ root  1   2   3   7   8   4   9  10   5  11  12   6  13  14 ]

   ro  [ C00                                                         ]
   01  [ C01 C11                                                     ]
   02  [ C01     C11                                                 ]
   03  [ C02 C12     C22                                             ]
   07  [ C03 C13     C23 C33                                         ]
   08  [ C03 C13     C23     C33                                     ]
   04  [ C02 C12                 C22                                 ]
   09  [ C03 C13                 C23 C33                             ]
   10  [ C03 C13                 C23     C33                         ]
   05  [ C02     C12                         C22                     ]
   11  [ C03     C13                         C23 C33                 ]
   12  [ C03     C13                         C23     C33             ]
   06  [ C02     C12                                     C22         ]
   13  [ C03     C13                                     C23 C33     ]
   14  [ C03     C13                                     C23     C33 ]

   The block entries can be easily deduced once the node-order is
   determined:

   - Block-Row nd has all entries from row pd(nd) from CORE-matrix C
   - every Block-Row has the corresponding diagonal entry
   - every Block-Row has also entries in columns corresponding to
     ancestor nodes of 'nd'. The column they are taken from in C
     corresponds to the period of the ancestor

   i.e. row corresponding to node 12: is period 3.
   Ancestors of 12 are root - 2 - 5 - 12
    => entries in columns corresponding to root, 2, 5, 12
       entries are C03, C13, C23, C33

   Any Algebra based on Schur-Complement is more efficient if the
   number of complicating columns is as small as possible.

   A column in any of the root periods that has entries only
   in blocks corresponding to the root-periods, can be treated
   separately as a further diagonal block:
   - determine columns in the root block, that do not have entries
     in blocks corresponding to diagonal elements (C0).
   - move these columns to the right of the root block.
   - find rows of the root block which are never referenced in
     the C0 columns; move them to the top (R0).
   =>
           C1 C0 D
       R0  B
       R1  B  A
       D   B  0  A
  */

  //
  // compute the dimensions of the submatrices
  //

  // nonzeros in RnkC, Diag parts for certain periods in core
  int *rnkc_nz_pd = new int[nPeriods];
  int *diag_nz_pd = new int[nPeriods];

  // count nonzeros in parts of big matrix
  for (i = 0; i < nPeriods; ++i) {

    rnkc_nz_pd[i] = diag_nz_pd[i] = 0;
    for (j = 0; j < cutoff; ++j)
      rnkc_nz_pd[i] += smps.getNzPeriod(i, j);
    for (j = cutoff; j < nPeriods; ++j)
      diag_nz_pd[i] += smps.getNzPeriod(i, j);

#ifdef DEBUG_GENERATE_SMPS
    printf("Nonzeros Diagon[%d]: %5d\tBorder[%d]: %5d\n",
           i, diag_nz_pd[i], i, rnkc_nz_pd[i]);
#endif
  }

  node = rootNode;

  do {

    perNode = node->level();
    blkNode = node->block();
    assert(blkNode < nBlocks + 1);

    for (i = 0; i < node->nLevels(); ++i) {
      rnkc_nz_blk[blkNode] += rnkc_nz_pd[perNode + i];
      diag_nz_blk[blkNode] += diag_nz_pd[perNode + i];
      diag_m_blk[blkNode]  += smps.getNRowsPeriod(perNode + i);
      diag_n_blk[blkNode]  += smps.getNColsPeriod(perNode + i);
    }
    /*
    printf("%d  %d  %d  %d ", node->name(), blkNode, perNode, node->nLevels());
    printf("* %6d %6d %6d %6d\n",
	   rnkc_nz_blk[blkNode], diag_nz_blk[blkNode],
	   diag_m_blk[blkNode],  diag_n_blk[blkNode]);
    */
  } while (node = node->next());

  rnkc_m_blk[0] = diag_m_blk[0];
  rnkc_n_blk[0] = diag_n_blk[0] - cldg0;
  diag_n_blk[0] = cldg0;

  // save dimensions for postprocess
  Ret.nRowsRnkc = diag_m_blk[0];
  Ret.nColsRnkc = rnkc_n_blk[0];
  Ret.nColsDiag = diag_n_blk[0];

  //
  // set the first row/col of each block (in big matrix before reordering)
  //
  int *f_rw_blk = new int[nBlocks + 2];
  int *f_cl_blk = new int[nBlocks + 2];

  f_rw_blk[0] = 0;
  f_cl_blk[0] = rnkc_n_blk[0];

  for (i = 1; i <= nBlocks; ++i) {
    rnkc_m_blk[i] = diag_m_blk[i];
    rnkc_n_blk[i] = rnkc_n_blk[0];
    f_rw_blk[i] = f_rw_blk[i - 1] + diag_m_blk[i - 1];
    f_cl_blk[i] = f_cl_blk[i - 1] + diag_n_blk[i - 1];
  }
  f_rw_blk[i] = f_rw_blk[i - 1] + diag_m_blk[i - 1];
  f_cl_blk[i] = f_cl_blk[i - 1] + diag_n_blk[i - 1];

  rnkc_nz_blk[0] -= nzdg0;
  diag_nz_blk[0] += nzdg0;

  //
  // allocate the sparse matrices for rank corrector and diagonal components
  //

  SparseSimpleMatrix *sparse;

  // diagonal algebras
  Algebra **Diagon = (Algebra**) malloc((nBlocks + 2) * sizeof(Algebra *));

  // rankcor part
  Algebra **Border = (Algebra**) malloc((nBlocks + 1) * sizeof(Algebra *));

  // bottom algebras (set to zero)
  Algebra **Bottom = (Algebra**) malloc((nBlocks + 1) * sizeof(Algebra *));

  // diagonal elements of Q
  Algebra **QDiag = (Algebra**) malloc((nBlocks + 2) * sizeof(Algebra *));

  char name[15];

#ifdef DEBUG_GENERATE_SMPS
#ifdef WITH_MPI
  if(IS_ROOT_PAR)
#endif
  printf("Dimensions and nonzeros of the deterministic equivalent blocks:\n"
         " blk | diagon (rows, cols)  nz | border (rows, cols)  nz ||"
	 " first row/col\n");
#endif /* DEBUG_GENERATE_SMPS */

  for (i = 0; i <= nBlocks; ++i) {

#ifdef DEBUG_GENERATE_SMPS
#ifdef WITH_MPI
  if(IS_ROOT_PAR)
#endif
    printf("%4d | (%6d, %6d) %6d | (%6d, %6d) %6d || %6d %6d\n", i,
           diag_m_blk[i], diag_n_blk[i], diag_nz_blk[i],
           rnkc_m_blk[i], rnkc_n_blk[i], rnkc_nz_blk[i],
           f_rw_blk[i],   f_cl_blk[i]);
#endif /* DEBUG_GENERATE_SMPS */

    // right-hand columns
    sprintf(name, "Border[%d]", i);
    sparse = NewSparseMatrix(rnkc_m_blk[i], rnkc_n_blk[i],
			     rnkc_nz_blk[i], name);
    sparse->cbf = (CallBackFunction) CallBackVoid;
    sparse->nb_el = 0;
    Border[i] = NewSparseSimpleAlgebra(sparse);

    // diagonal entries
    sprintf(name, "Diagon[%d]", i);
    sparse = NewSparseMatrix(diag_m_blk[i], diag_n_blk[i],
			     diag_nz_blk[i], name);
    Diagon[i] = NewSparseSimpleAlgebra(sparse);
    sparse->cbf = (CallBackFunction) CallBackVoid;
    sparse->nb_el  = 0;
    sparse->nb_col = 0;

    // bottom part
    sprintf(name, "Bottom[%d]", i);
    sparse = NewSparseMatrix(0, diag_n_blk[i], 0, name);
    sparse->cbf = (CallBackFunction) CallBackVoid;
    Bottom[i] = NewSparseSimpleAlgebra(sparse);

    // Q diagonal part
    sprintf(name, "Q_Diag[%d]", i);
    sparse = NewSparseMatrix(diag_n_blk[i], diag_n_blk[i], 0, name);
    sparse->cbf = (CallBackFunction) CallBackVoid;
    QDiag[i] = NewSparseSimpleAlgebra(sparse);
  }

  // last diagonal corresponding to the bottom blocks
  sparse = NewSparseMatrix(0, rnkc_n_blk[0], 0, "DiagPart");
  sparse->cbf = (CallBackFunction) CallBackVoid;
  Diagon[nBlocks + 1] = NewSparseSimpleAlgebra(sparse);

  sparse = NewSparseMatrix(rnkc_n_blk[0], rnkc_n_blk[0], 0, "QDiagPart");
  sparse->cbf = (CallBackFunction) CallBackVoid;
  QDiag[nBlocks + 1] = NewSparseSimpleAlgebra(sparse);

  //
  // build the deterministic equivalent column by column
  //

  node = rootNode;
  blkNode = 0;

  // counter of columns added to Rnkc so far
  int ncol_rc = 0;

  // first column in current block
  int fColBlk = 0;

  int cIndex = 0, idxCol;
  Algebra **Array;

  // for all columns in the deterministic equivalent
  for (int col = 0; col < ttn; ++col) {

    int lastCol = 0;
    perNode = node->level();

    for (j = 0; j < node->nLevels(); ++j)
      lastCol += smps.getNColsPeriod(perNode + j);

    if (col - fColBlk >= lastCol) {

      // first col of current node in big matrix
      fColBlk = col;

      // current block
      node = node->next();
      assert(node != NULL);

      // period of current node
      perNode = node->level();

      // current part of big matrix
      blkNode = node->block();

      // index in the core matrix
      cIndex = data.clpnts[smps.getBegPeriodCol(perNode)];

#ifdef DEBUG_GENERATE_SMPS
      printf("Column %4d - node now %3d  stage now %d  block now %d\n",
             col, node->name(), perNode, blkNode);
#endif
    }

    // corresponding column in the core matrix
    int coreCol = col - fColBlk + smps.getBegPeriodCol(perNode);
    assert(coreCol <= smps.getCols());

    //
    // initialise the column in the current block
    //

    // border column
    if (blkNode == 0 && is_col_diag[coreCol] == 0) {

      for (j = 0; j <= nBlocks; ++j) {
        sparse = (SparseSimpleMatrix *) Border[j]->Matrix;
        sparse->col_beg[ncol_rc] = sparse->nb_el;
      }
      sparse = (SparseSimpleMatrix *) Border[0]->Matrix;
      idxCol = ncol_rc;
      ++ncol_rc;
      Array = Border;
    }

    // diagonal column
    else {
      sparse = (SparseSimpleMatrix *) Diagon[blkNode]->Matrix;
      sparse->col_beg[sparse->nb_col] = sparse->nb_el;
      idxCol = sparse->nb_col;
      sparse->nb_col++;
      Array = Diagon;
    }

    int offset = node->firstRow() - smps.getBegPeriodRow(perNode)
      - f_rw_blk[blkNode];

    // number of nonzeros in the current column of the core matrix
    int nnzCol = data.clpnts[coreCol + 1] - data.clpnts[coreCol];

    while (nnzCol > 0) {

      int row = data.rwnmbs[cIndex];

      // stop early before starting a new period
      if (row >= smps.getBegPeriodRow(perNode + node->nLevels()))
	break;

      // objective row
      if (row == smps.getObjRowIndex()) {
	// nothing to do
      }

      else {
        assert(sparse->nb_el < sparse->max_nb_el);
        sparse->element[sparse->nb_el] = data.acoeff[cIndex];
        sparse->row_nbs[sparse->nb_el] = row + offset;

#ifdef DEBUG_GENERATE_SMPS
	printf(" node %2d :> %2d - %2d + %2d - %2d  ", node->name(),
	       row, smps.getBegPeriodRow(perNode),
	       node->firstRow(), f_rw_blk[blkNode]);
	printf(":: % lf  (row %d)\n", sparse->element[sparse->nb_el],
	       sparse->row_nbs[sparse->nb_el]);
#endif

        assert(sparse->row_nbs[sparse->nb_el] >= 0);
        assert(sparse->row_nbs[sparse->nb_el] < sparse->nb_row);

	sparse->nb_el++;
	sparse->col_len[idxCol]++;
      }

      ++cIndex;
      --nnzCol;
    }

    // write all parts of this column in the correct border and diagonal
    // matrices by scanning through the node and its children
    copyLinkingBlocks(smps, Array, data, node, f_rw_blk,
                      idxCol, cIndex, nnzCol);
  }

  assert(Ret.nColsRnkc == ncol_rc);

  // write the final pointer for the last+1 column
  for (j = 0; j <= nBlocks; ++j) {
    sparse = (SparseSimpleMatrix *) Border[j]->Matrix;
    sparse->col_beg[Ret.nColsRnkc] = sparse->nb_el;
    sparse = (SparseSimpleMatrix *) Diagon[j]->Matrix;
    sparse->col_beg[sparse->nb_col] = sparse->nb_el;
  }

  Ret.b = NewDenseVector(ttm, "RHS");
  Ret.c = NewDenseVector(ttn, "Obj");
  Ret.l = NewDenseVector(ttn, "LB");
  Ret.u = NewDenseVector(ttn, "UB");

  // setup the right-hand side
  setupRhs(smps, &Ret);

  // setup objective and bounds
  setupObjective(smps, &Ret);

  // apply scenario corrections
  applyScenarios(tree, &Ret, Diagon, Border, f_rw_blk, f_cl_blk);

  // reorder objective, bounds and column names
  reorderObjective(tree, &Ret, rnkc_n_blk[0]);

  // set up the deterministic equivalent as a DblBordDiagAlgebra
  DblBordDiagSimpleMatrix *MA =
    NewDblBordDiagSimpleMatrix(nBlocks + 1, Bottom, Border, Diagon, "SmpsA");
  BlockDiagSimpleMatrix *MQ =
    NewBlockDiagSimpleMatrix(nBlocks + 2, QDiag, "Q_main");

  Ret.AlgA = NewDblBordDiagSimpleAlgebra(MA);
  Ret.AlgQ = NewBlockDiagSimpleAlgebra(MQ);

  // clean up
  delete[] rnkc_m_blk;
  delete[] rnkc_n_blk;
  delete[] rnkc_nz_blk;
  delete[] rnkc_nz_pd;
  delete[] diag_m_blk;
  delete[] diag_n_blk;
  delete[] diag_nz_blk;
  delete[] diag_nz_pd;
  delete[] f_rw_blk;
  delete[] f_cl_blk;

  return 0;
}

/**
 *  Copy the elements in the linking blocks.
 *
 *  Go through the given node and its children (recursively): for each node
 *  append the relevant period part of the core matrix column to the current
 *  column.
 *
 *  @param smps:
 *         The smps instance
 *  @param Array:
 *         Pointer to an array of algebras
 *  @param data:
 *         The sparse representation of the core matrix
 *  @param node:
 *         Node in the deterministic equivalent at which to start
 *  @param f_rw_blk:
 *         First row of each block in the deterministic equivalent
 *  @param idxCol:
 *         Index of the current column in the deterministic equivalent
 *  @param cIndex:
 *         Index of the sparse core data
 *  @param nnzCol:
 *         Remaining number of nonzeros in the current column
 */
int copyLinkingBlocks(Smps &smps, Algebra **Array, const SparseData &data,
                      const Node *node, int *f_rw_blk,
                      int &idxCol, int &cIndex, int &nnzCol) {

  // store the current row index of core
  const int sIndex = cIndex;

  // store the number of remaining nonzeros in the column
  const int snzCol = nnzCol;

  // the offset in row numbers between core and deterministic equivalent
  int offset;

  // the period of the node
  const int period = node->level();

  SparseSimpleMatrix *sparse;

  // copy the linking blocks from the current period
  for (int chd = 0; chd < node->nChildren(); ++chd) {

    const Node *child = node->getChild(chd);
    sparse = (SparseSimpleMatrix *) Array[child->block()]->Matrix;

    // determine the offset in row numbers for the children nodes
    offset = child->firstRow() - smps.getBegPeriodRow(child->level()) -
      f_rw_blk[child->block()];

    // restore the row index of core and the number of nonzeros
    cIndex = sIndex;
    nnzCol = snzCol;

    // copy all the remaining nonzeros in the current column
    while (nnzCol > 0) {

      // an element in a linking block may belong to the next period,
      // in which case it has to be copied as many times as there are
      // children of the current node, or it belongs to a period after
      // the next (not staircase structure), in which case it has to be
      // copied as many times as there are nodes in that period.
      // this second case is accomplished by recursion on the children
      // of the current node

      int row = data.rwnmbs[cIndex];

      // this element belongs to the next period
      if (row < smps.getBegPeriodRow(period + 2)) {
	assert(sparse->nb_el < sparse->max_nb_el);
        sparse->element[sparse->nb_el] = data.acoeff[cIndex];
	sparse->row_nbs[sparse->nb_el] = row + offset;

#ifdef DEBUG_GENERATE_SMPS
        printf(" node %2d :> %2d - %2d + %2d - %2d  ", child->name(),
               row, smps.getBegPeriodRow(child->level()),
               child->firstRow(), f_rw_blk[child->block()]);
        printf(":: % lf  (row %d)\n", sparse->element[sparse->nb_el],
               sparse->row_nbs[sparse->nb_el]);
#endif

	assert(sparse->row_nbs[sparse->nb_el] >= 0);
	assert(sparse->row_nbs[sparse->nb_el] < sparse->nb_row);

	sparse->nb_el++;
	sparse->col_len[idxCol]++;

	++cIndex;
	--nnzCol;
      }

      // this element belongs to a period after the next
      else
	copyLinkingBlocks(smps, Array, data, child, f_rw_blk,
                          idxCol, cIndex, nnzCol);
    }
  }

  return 0;
}

/** Set up the right-hand side */
void setupRhs(const Smps &smps, SmpsReturn *Ret) {

  int firstRowNode, begRowPeriod;
  DenseVector *rhs = Ret->b;
  const Node *node = Ret->rootNode;

  // leave immediately if there is no root node
  if (!node)
    return;

  // generate the row names for the deterministic equivalent
  Ret->rownames = smps.getRowNames();

  // for all nodes in the tree in order
  do {

    firstRowNode = node->firstRow();
    begRowPeriod = smps.getBegPeriodRow(node->level());

    // copy the information for this node
    for (int i = 0; i < node->nRows(); ++i) {

      rhs->elts[firstRowNode + i] = smps.getRhs(begRowPeriod + i);
    }

  } while (node = node->next());
}

/** Set up the objective and the bounds */
void setupObjective(const Smps &smps, SmpsReturn *Ret) {

  int firstColNode, begColPeriod;
  DenseVector *obj = Ret->c, *lob = Ret->l, *upb = Ret->u;
  const Node *node = Ret->rootNode;

  // leave immediately if there is no root node
  if (!node)
    return;

  // generate the column names for the deterministic equivalent
  Ret->colnames = smps.getColNames();

  // copy the objective row from the core matrix
  double *coreObj = smps.getObjRow();

  // for all nodes in the tree
  do {

    firstColNode = node->firstCol();
    begColPeriod = smps.getBegPeriodCol(node->level());
    double probNode = node->probNode();

    // copy the information for this node
    for (int i = 0; i < node->nCols(); ++i) {

      // copy the objective coefficients weighted by probability of the node
      obj->elts[firstColNode + i] = probNode * coreObj[begColPeriod + i];

      // copy the bounds
      lob->elts[firstColNode + i] = smps.getLowerBound(begColPeriod + i);
      upb->elts[firstColNode + i] = smps.getUpperBound(begColPeriod + i);
    }

  } while (node = node->next());

  // clean up
  delete[] coreObj;
}

/**
 *  Apply the scenarios changes.
 *
 *  Every scenario is a path from the root node to one of the leaf
 *  nodes. The scenario change data has for all scenarios a list of
 *  changes that need to be applied. These changes can apply to any
 *  period of the tree (not just the leaf nodes).
 *
 *  @param tree
 *         The tree for which we are building the deterministic equivalent
 *  @param Ret
 *         The SmpsReturn structure of the problem
 *  @param Diagon
 *         Array of diagonal algebras
 *  @param Border
 *         Array of right column algebras already in final ordering: rows
 *         are ordered as [D1, ..., Dn, D0], columns as [RnkCor, D0]
 *  @param f_rw_blk
 *         First row of each block in the deterministic equivalent
 *  @param f_cl_blk
 *         First column of each block in the deterministic equivalent
 *
 *  @note
 *  The entryCol[] and entryRow[] vectors are in FORTRAN numbering.
 *
 *  Need: for each node, start of row/col within its DetEquivMatrix blck
 *   -   for diagonal blocks can take row/col from f_cl_nd, f_rw_nd and
 *       substract starting row/col for seed node of this block
 *   -   for nondiag node need to do breadth first count
 *       the correct column within the RnckD0 part is obtained by looping
 *       through all columns within RnckD0 (breadth first over nodes)
 *       and counting how many end up in Rnck and D0
 */
int SmpsOops::applyScenarios(const SmpsTree &tree, SmpsReturn *Ret,
                             Algebra **Diagon, Algebra **Border,
			     int *f_rw_blk, int *f_cl_blk) {

  int period, lastPd, firstEntry, lastEntry;
  int row, col, pdr, pdc;
  int k, iblk, jblk;
  bool found;
  const int *entryRow = smps.getEntryRow();
  const int *entryCol = smps.getEntryCol();
  const double *entryVal = smps.getEntryVal();
  const Node *node = Ret->rootNode, *scNode;
  SparseSimpleMatrix *sparse;

  // We start from the leaf nodes and traverse the tree up to the root,
  // applying the changes corresponding to the nodes in the scenario
  // only if the correction refers to the same period of the node.

  do {

    // skip the non leaf nodes
    if (node->nChildren() > 0)
      continue;

    scNode = node;
    period = node->level();

    // traverse all nodes from this leaf up to the root
    while (period >= 0) {

      // last period covered by this node
      lastPd = period + scNode->nLevels();

      // scenario the current node belongs to
      int scen = scNode->scenario();

      // index of the first change
      firstEntry = smps.getFirstEntryScen(scen) - 1;
      lastEntry  = firstEntry + smps.getLengthScen(scen);

#ifdef DEBUG_SCEN
      printf("Node %d: scen %d (entries from %d to %d), period %d to %d\n",
	     scNode->name(), scen + 1, firstEntry, lastEntry - 1,
	     period, lastPd);
#endif

      // for all changes affecting this scenario
      for (int corr = firstEntry; corr < lastEntry; ++corr) {

	// row and column of core affected by the change
	pdr = smps.getRowPeriod(entryRow[corr] - 1);
	pdc = smps.getColPeriod(entryCol[corr] - 1);

	// if the change affects the objective
	if (pdr < 0 && pdc >= period && pdc < lastPd) {

	  col = scNode->firstCol() + entryCol[corr] - 1
	    - smps.getBegPeriodCol(period);
	  assert(col <= tree.getTotCols());

	  Ret->c->elts[col] = entryVal[corr] * scNode->probNode();

#ifdef DEBUG_SCEN
	  printf("   Obj entry: core col %d, det.eq. col %d, (%f)\n",
		 entryCol[corr], col, entryVal[corr]);
#endif
	}

	// if the change affects the rhs
	else if (pdc < 0 && pdr >= period && pdr < lastPd) {

	  row = scNode->firstRow() + entryRow[corr] - 1
	    - smps.getBegPeriodRow(period);
	  assert(row <= tree.getTotRows());

	  Ret->b->elts[row] = entryVal[corr];

#ifdef DEBUG_SCEN
	  printf("   Rhs entry: core row %d, det.eq. row %d, (%f)\n",
		 entryRow[corr], row, entryVal[corr]);
#endif
	}

	// if the change affects the matrix
	else if (pdc >= 0 && pdr >= 0 && pdr == period) {

	  // we cannot deal with above-diagonal entries
	  assert((pdr == pdc) || (pdr == pdc + 1));

	  // indices in the deterministic equivalent
	  row = scNode->firstRow() + entryRow[corr] - 1
	    - smps.getBegPeriodRow(pdr);
	  col = scNode->firstCol() + entryCol[corr] - 1
	    - smps.getBegPeriodCol(pdc);

	  // adjust the column index if the change is in a column that belongs
	  // to the previous period
	  if (pdr != pdc)
	    col += scNode->parent()->firstCol() - scNode->firstCol();

	  assert(row <= tree.getTotRows());
	  assert(col <= tree.getTotCols());

#ifdef DEBUG_SCEN
	  printf("   Mtx entry: core row %3d col %3d, det.eq. row %3d col %3d",
		 entryRow[corr], entryCol[corr], row, col);
#endif

	  // work out the row and column block of the correction
	  int rowBlock = 0;
	  while (row >= f_rw_blk[rowBlock + 1]) ++rowBlock;
	  jblk = row - f_rw_blk[rowBlock];
	  int colBlock = 0;
	  while (col >= f_cl_blk[colBlock + 1]) ++colBlock;
	  iblk = col - f_cl_blk[colBlock];

	  assert(rowBlock <= nBlocks);
	  assert(colBlock <= nBlocks);
	  found = false;

	  // if change is in the RankCor part of the deterministic equivalent
	  if (colBlock == 0) {

	    // work out if this column is in the RankCor part or not by
	    // looping through all RankCor columns and counting

	    // columns in rnk and in d0 so far
	    int iblkrnc = 0, iblkd0 = 0;

	    // node of the current column and corresponding period
	    const Node *nd = Ret->rootNode;
	    int ndPd = nd->level();
	    int bigCol, coreCol;

	    // find the column in the core matrix that this belongs to
	    for (bigCol = 0, coreCol = 0; bigCol < col; ++bigCol, ++coreCol) {

	      int perNode = nd->level();

	      // if this is past the last column in this node
	      if (coreCol >= smps.getBegPeriodCol(perNode)) {
		nd = nd->next();
		ndPd = perNode;
		coreCol = smps.getBegPeriodCol(ndPd);
	      }
	      if (Ret->is_col_diag[coreCol])
		++iblkd0;
	      else
		++iblkrnc;
	    }

	    if (Ret->is_col_diag[coreCol]) {
	      sparse = (SparseSimpleMatrix *) Diagon[0]->Matrix;
	      iblk = iblkd0;
	    } else {
	      sparse = (SparseSimpleMatrix *) Border[rowBlock]->Matrix;
	      iblk = iblkrnc;
	    }
	  }

	  // if the change is in the diagonal part
	  else {

	    // ensure we are in a diagonal block
	    assert(colBlock == rowBlock);

	    sparse = (SparseSimpleMatrix *) Diagon[colBlock]->Matrix;
	  }

	  // apply the change
	  for (k = sparse->col_beg[iblk]; k < sparse->col_beg[iblk + 1]; ++k) {
	    if (sparse->row_nbs[k] == jblk) {
#ifdef DEBUG_SCEN
	      printf(" (%g -> %g)\n", sparse->element[k], entryVal[corr]);
#endif
	      sparse->element[k] = entryVal[corr];
	      found = true;
	    }
	  }
	  if (!found) {
	    printf("\nEntry not found\n");
	    exit(1);
	  }
	}
#ifdef DEBUG_SCEN
	// the change is in a different period, so we don't apply it now
	else {
	  int pchange = (pdr > pdc) ? pdr: pdc;
	  printf("   Ignored change in period %d (now is period %d)\n",
		 pchange + 1, period + 1);
	  assert(pchange != period);
	}
#endif
      }

      // walk up the tree to the parent node at the previous period
      scNode = scNode->parent();
      period--;
    }

  } while (node = node->next());

  return 0;
}

/**
 *  Reorder objective and bounds.
 *
 *  Columns have been set up in the order:
 *     [Block-0 Diag-1 ... Diag-n]
 *  However, in the deterministic equivalent Block-0 is split in RnkC + Diag-0
 *  (some columns from the first period are pushed into the diagonal part).
 *  Here columns are changed to match the ordering in the deterministic
 *  equivalent:
 *     [Diag-0 Diag-1 ... Diag-n RnkC]
 */
void SmpsOops::reorderObjective(const SmpsTree &tree, SmpsReturn *Ret,
				const int rnkn) {

  int col, coreCol, firstColDiag, firstColNode;
  int nb_el = 0;
  int ttn = tree.getTotCols();

  const Node *node = Ret->rootNode;
  if (!node)
    return;

  DenseVector *obj = Ret->c, *lob = Ret->l, *upb = Ret->u;
  int *is_col_diag = Ret->is_col_diag;
  char **colnames  = Ret->colnames;

  double *objCopy = (double *) malloc(ttn * sizeof(double));
  double *upbCopy = (double *) malloc(ttn * sizeof(double));
  double *lobCopy = (double *) malloc(ttn * sizeof(double));
  char  **clnCopy = NULL;
  if (colnames)
    clnCopy = (char **) malloc(ttn * sizeof(char *));

  // find the first node that will go in a diagonal block
  while (node->level() < cutoff)
    node = node->next();

  // set the first column in diagonal block
  firstColDiag = node->firstCol();

  // copy objective and bounds into temporary arrays
  for (int i = 0; i < ttn; ++i) {
    objCopy[i] = obj->elts[i];
    lobCopy[i] = lob->elts[i];
    upbCopy[i] = upb->elts[i];
    if (colnames)
      clnCopy[i] = colnames[i];
  }

  // copy the diagonal elements from the first block (Diag-0)
  firstColNode = 0;
  for (col = 0, node = Ret->rootNode; col < firstColDiag; ++col) {

    // check if the current column belongs to the next node
    if (col - firstColNode >= node->nCols()) {
      firstColNode = col;
      node = node->next();
    }

    // find the corresponding column in the core matrix
    coreCol = col - firstColNode + smps.getBegPeriodCol(node->level());
    assert(coreCol <= smps.getCols());
    if (is_col_diag[coreCol] == 1) {
      obj->elts[nb_el] = objCopy[col];
      lob->elts[nb_el] = lobCopy[col];
      upb->elts[nb_el] = upbCopy[col];
      if (colnames)
	colnames[nb_el] = clnCopy[col];
      ++nb_el;
    }
  }

  assert(nb_el == firstColDiag - rnkn);

  // copy the other diagonal entries  (Diag-1 ... Diag-n)
  for (col = firstColDiag; col < ttn; ++col) {
    obj->elts[nb_el] = objCopy[col];
    lob->elts[nb_el] = lobCopy[col];
    upb->elts[nb_el] = upbCopy[col];
    if (colnames)
      colnames[nb_el] = clnCopy[col];
    ++nb_el;
  }

  // copy the RankCorrector entries from the first block
  firstColNode = 0;
  for (col = 0, node = Ret->rootNode; col < firstColDiag; ++col) {

    // check if the current column belongs to the next node
    if (col - firstColNode >= node->nCols()) {
      firstColNode = col;
      node = node->next();
    }

    // find the corresponding column in the core matrix
    coreCol = col - firstColNode + smps.getBegPeriodCol(node->level());
    assert(coreCol <= smps.getCols());
    if (is_col_diag[coreCol] == 0) {
      obj->elts[nb_el] = objCopy[col];
      lob->elts[nb_el] = lobCopy[col];
      upb->elts[nb_el] = upbCopy[col];
      if (colnames)
	colnames[nb_el] = clnCopy[col];
      ++nb_el;
    }
  }

  assert(nb_el == ttn);

  // clean up
  free(objCopy);
  free(lobCopy);
  free(upbCopy);
  free(clnCopy);
}

/**
 *  Copy a Vector into a breadth-first ordered DenseVector.
 *
 *  Copies a Vector as used by OOPS into a DenseVector corresponding to a
 *  breadth-first ordering of the scenario tree.
 */
void SmpsOops::VectorToSmpsDense(Vector *x, DenseVector *dx,
                                 const SmpsReturn &Ret, const int rowcol) {

  SetExactVector(x);
  CopyToDenseVector(x, dx);

  // reorder according to the SMPS breadth-first order
  if (rowcol == ORDER_COL)
    backOrderColVector(smps, Ret, dx->elts);
  else
    backOrderRowVector(Ret, dx->elts);
}

/**
 *  Copy a breadth-first ordered DenseVector into a Vector.
 *
 *  Copies a DenseVector corresponding to a breadth-first ordering of the
 *  scenario tree into a Vector in the mixed ordering used by OOPS.
 *
 *  Reorders the elements according to the reordering used by OOPS:
 *  first period entries are placed at the end rather than at the beginning,
 *  those columns that are not linking periods are placed in separate
 *  diagonal block, rather than in the RankCor block.
 */
void SmpsOops::SmpsDenseToVector(DenseVector *dx, Vector *x,
                                 const SmpsReturn &Ret, const int rowcol) {

  // should attempt to leave the original element order intact

  // go for memory saving option
  if (rowcol == ORDER_COL)
    forwOrderColVector(smps, Ret, dx->elts);
  else
    forwOrderRowVector(Ret, dx->elts);

  CopyDenseToVector(dx, x);

  // and reverse the order, to leave the dense vector intact
  if (rowcol == ORDER_COL)
    backOrderColVector(smps, Ret, dx->elts);
  else
    backOrderRowVector(Ret, dx->elts);
}

/**
 *  The column vectors are setup in the order
 *   [D0 D1 D2 ... Dn-1 Rnk]
 *  ([D0 Rnk] is the RankCorrector and D1,...,Dn-1 are the 'proper' diagonal
 *  blocks).
 *
 *  This routine re-orders them into the SMPS order:
 *   [(D0Rnk) D1 D2 ... Dn-1]
 *  where (D0Rnk) is the reordered [D0 Rnk] using the original
 *  order of these columns in the core matrix.
 *
 *  Needed from main method (passed through SmpsReturn *Ret)
 *  - nColsRnkc: number of columns in actual rankcor (Rnk) part
 *  - nColsDiag: number of columns in diag rankcor (D0) part
 */
void backOrderColVector(const Smps &smps, const SmpsReturn &Ret, double *x) {

  const Node *node = Ret.rootNode;

  // total number of columns in RankCor (D0|Rnk)
  const int ncol_ttrc = Ret.nColsRnkc + Ret.nColsDiag;
  const int ttn = Ret.c->dim;

  double *dtmp  = new double[ttn];

  for (int i = 0; i < ttn; ++i)
    dtmp[i] = x[i];

  // copy entries into the combined rankcor slot

  // set the start of the D0 and Rnk blocks
  int nx_col_d0 = 0;                      // next column to take from D0
  int nx_col_rc = ttn - Ret.nColsRnkc;    // next column to take from Rnk
  int currPer   = 0;                      // period of the current column

#ifdef DEBUG_SMPS_ORDER
  printf("BACKORDCOL: node %3d (per %d): core cols %d--%d\n",
	 node->name(), currPer,
	 0, smps.getBegPeriodCol(currPer + 1) - 1);
#endif

  // start writing the (D0Rnk) part
  for (int col = 0, coreCol = 0; col < ncol_ttrc; ++col, ++coreCol) {

    // find the column in the original core that this belongs to
    if (coreCol >= smps.getBegPeriodCol(currPer + 1)) {

      // if this is past the last column in this node
      node = node->next();
      currPer = node->level();
      coreCol = smps.getBegPeriodCol(currPer);

#ifdef DEBUG_SMPS_ORDER
      printf("BACKORDCOL: node %3d (per %d): core cols %d--%d\n",
      printf("BACKORDCOL: node %3d: pd=%2d f_cl_core=%d l_cl_core=%d\n",
	     node->name(), currPer,
	     coreCol, smps.getBegPeriodCol(currPer + 1) - 1);
#endif
    }

    // copy the value depending on whether it is in D0 or Rnk
    if (Ret.is_col_diag[coreCol])
      x[col] = dtmp[nx_col_d0++];
    else
      x[col] = dtmp[nx_col_rc++];
  }

  // copy the remaining columns in order, beginning after rankcor
  int offset = ncol_ttrc;

  // nx_col_d0 points to the next col that should be copied from main part
  for (int col = 0; col < ttn - ncol_ttrc; ++col) {
    x[col + offset] = dtmp[nx_col_d0++];
  }

  // clean up
  delete[] dtmp;
}

/**
 *  The column vectors are setup in OOPS in the order
 *   [D0 D1 D2 ... Dn-1 Rnk]
 *  ([D0 Rnk] is the RankCorrector and D1,...,Dn-1 are the 'proper' diagonal
 *  blocks).
 *
 *  This routine re-orders them from the SMPS order:
 *   [(D0Rnk) D1 D2 ... Dn-1]
 *
 *  In both cases (D0Rnk) is the reordered [D0 Rnk] using the original
 *  order of these columns in the core matrix.
 *
 *  Needed from main method (passed through SmpsReturn *Ret)
 *  - nColsRnkc: number of columns in actual rankcor (Rnk) part
 *  - nColsDiag: number of columns in diag rankcor (D0) part
 */
void forwOrderColVector(const Smps &smps, const SmpsReturn &Ret, double *x) {

  const Node *node = Ret.rootNode;

  // total number of columns in RankCor (D0|Rnk)
  const int ncol_ttrc = Ret.nColsRnkc + Ret.nColsDiag;
  const int ttn = Ret.c->dim;
  double *dtmp  = new double[ttn];

  for (int i = 0; i < ttn; ++i)
    dtmp[i] = x[i];

  // copy entries from the combined rankcor slot into OOPS RankCor and
  // first diagonal

  // set the start of the D0 and Rnk blocks
  int nx_col_d0 = 0;                      // next column in D0
  int nx_col_rc = ttn - Ret.nColsRnkc;    // next column in Rnk
  int currPer   = 0;                      // period of the current column

#ifdef DEBUG_SMPS_ORDER
  printf("FORWORDCOL: node %3d (per %d): core cols %d--%d\n",
	 node->name(), currPer,
	 0, smps.getBegPeriodCol(currPer + 1) - 1);
#endif

  // loop through all columns in the D0Rnk block of original vector
  for (int col = 0, coreCol = 0; col < ncol_ttrc; ++col, ++coreCol) {

    // find the column in the original core that this belongs to
    if (coreCol >= smps.getBegPeriodCol(currPer + 1)) {

      // if this is past the last column in this node
      node = node->next();
      currPer = node->level();
      coreCol = smps.getBegPeriodCol(currPer);

#ifdef DEBUG_SMPS_ORDER
      printf("FORWORDCOL: node %3d (per %d): core cols %d--%d\n",
	     node->name(), currPer,
	     coreCol, smps.getBegPeriodCol(currPer + 1) - 1);
#endif
    }

    // copy the value depending on whether it is in D0 or Rnk
    if (Ret.is_col_diag[coreCol])
      x[nx_col_d0++] = dtmp[col];
    else
      x[nx_col_rc++] = dtmp[col];
  }

  // copy the remaining columns in order, beginning after rankcor
  int offset = ncol_ttrc;

  // nx_col_d0 points to the next entry that should be copied into
  for (int col = 0; col < ttn - ncol_ttrc; ++col) {
    x[nx_col_d0++] = dtmp[col + offset];
  }

  // clean up
  delete[] dtmp;
}

/**
 *  Reorder the rows into SMPS breadth-first order.
 *
 *  @param Ret:
 *         Information about the problem with respect to which the
 *         reordering should be done
 *  @param x:
 *         The vector to be reordered
 */
void backOrderRowVector(const SmpsReturn &Ret, double *x) {

  int i;
  const int ttm = Ret.b->dim, nRowsRnkc = Ret.nRowsRnkc;
  double *dtmp = new double[ttm];

  for (i = 0; i < ttm; ++i)
    dtmp[i] = x[i];

  // first copy the rows from the rankcor part (at end of matrix)
  for (i = 0; i < nRowsRnkc; ++i) {
    x[i] = dtmp[ttm - nRowsRnkc + i];
  }

  // now copy the rest
  for (i = nRowsRnkc; i < ttm; ++i) {
    x[i] = dtmp[i - nRowsRnkc];
  }

  // clean up
  delete[] dtmp;
}

/**
 *  Reorder the rows into the order used by OOPS.
 *
 *  Order the rows from the SMPS breadth-first order into the order used
 *  internally by OOPS (the only difference is that OOPS has the rows
 *  correspoding to the first 'cutoff' periods at the end, whereas the SMPS
 *  breadth-first order has them at the beginning).
 *
 *  @param Ret:
 *         Information about the problem with respect to which the
 *         reordering should be done
 *  @param x:
 *         The vector to be reordered
 */
void forwOrderRowVector(const SmpsReturn &Ret, double *x) {

  int i;
  const int ttm = Ret.b->dim, nRowsRnkc = Ret.nRowsRnkc;
  double *dtmp = new double[ttm];

  for (i = 0; i < ttm; ++i)
    dtmp[i] = x[i];

  // the rankcor rows are at the beginning of the vector, they should
  // be re-ordered to the end
  for (i = 0; i < nRowsRnkc; ++i) {
    x[ttm - nRowsRnkc + i] = dtmp[i];
  }

  // now copy the rest
  for (i = nRowsRnkc; i < ttm; ++i) {
    x[i - nRowsRnkc] = dtmp[i];
  }

  // clean up
  delete[] dtmp;
}
