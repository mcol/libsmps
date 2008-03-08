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

#include <assert.h>
#include "SmpsOops.h"


static void
setupRhs(const Smps &smps, SmpsReturn *Ret);

static void
setupObjective(const Smps &smps, SmpsReturn *Ret);


/**
 *  Generate the OOPS structures for the SMPS problem.
 *
 *  Uses the information from the SMPS files create the deterministic
 *  equivalent problem.
 *
 *  All nodes in periods < level are treated in a RankCorrector, nodes
 *  in period == level are seeds for diagonal blocks.
 *
 *  Some extra analysis is done to:
 *  - remove columns from RankCor that do not affect any periods >= level,
 *    these can be (and are) treated in an additional diagonal block (D-0)
 *    FIXME: should also separate rows from the first block of RankCor that
 *           do not affect columns in D-0
 *
 *  The subroutine procedes as follows:
 *
 *  (1a) Analyse which columns of the first block (periods < level) can be
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
SmpsReturn* SmpsOops::generateSmps(const Node *rootNode) {

  // ensure that the root node exists
  if (!rootNode)
    return NULL;

  // whole matrices A/Q
  Algebra *AlgA, *AlgQ;

  SmpsReturn *Ret = new SmpsReturn;

  int i, j, k;
  const Node *node;

  // store the root node
  Ret->rootNode = rootNode;

  // col and nonzeros in new D-0 matrix
  int nzdg0, cldg0;
  int *is_col_diag;   /* indicates columns in first periods that should be
			 in RankCor part but are not used outside these
			 periods. Hence they are taken out of the RankCor
			 and form an additional diagonal block */

  // node/period/block of current row/col
  int cu_nd_cl, cu_pd_rw, cu_pd_cl, cu_blk_cl;

  // counter of cols added to RnkC so far
  int ncol_rc;

  // first col in current block
  int b_cu_blk_cl;

  // number of periods
  const int nPeriods = smps.getPeriods();

  // pointer to start of period information in core for this column
  int *p_pd_rw = new int[nPeriods + 1];

  // dimensions and nonzeros of rank corrector blocks
  int *rnkc_m_blk  = new int[nBlocks + 1];
  int *rnkc_n_blk  = new int[nBlocks + 1];
  int *rnkc_nz_blk = new int[nBlocks + 1];

  // dimensions and nonzeros of diagonal blocks
  int *diag_m_blk  = new int[nBlocks + 1];
  int *diag_n_blk  = new int[nBlocks + 1];
  int *diag_nz_blk = new int[nBlocks + 1];

  // initialise to zero
  memset(diag_m_blk,  0, (nBlocks + 1) * sizeof(int));
  memset(diag_n_blk,  0, (nBlocks + 1) * sizeof(int));
  memset(diag_nz_blk, 0, (nBlocks + 1) * sizeof(int));
  memset(rnkc_nz_blk, 0, (nBlocks + 1) * sizeof(int));

  // nonzeros in RnkC, Diag parts for certain periods in core
  int *rnkc_nz_pd = new int[nPeriods];
  int *diag_nz_pd = new int[nPeriods];

  // dimensions of the deterministic equivalent
  int ttm = smps.getTotRows();
  int ttn = smps.getTotCols();

  // first row/col in the diagonal part of the deterministic equivalent
  int f_rwdiag = smps.getBegPeriodRow(level);
  int f_cldiag = smps.getBegPeriodCol(level);

  int objRow = smps.getObjRowIndex();

  const SparseData data = smps.getSparseData();

  SparseSimpleMatrix *sparse;

  printf(" --------------- generateSmps --------------\n");

  /*
     Scan through the columns of RankCorrector block to see which can go in
     the diagonal.
     For all row periods <= level, count the number of nonzeros in columns
     that can go in the diagonal part (used to obtain the number of nonzeros
     in diagonal entry).
     For all col periods <= level, count the number of columns that
     can go in the diagonal part (used to obtain col of diag).

     Columns within periods <= level are designated to go into the Border
     part of the matrix. However some of them have no entries in rows
     corresponding to other time periods, and they can be treated differently.

     This code identifies these columns (is_col_diag[i]) and counts them.
     Since periods other than the first one might be in the rank corrector,
     and non-first period columns of the core matrix will be repeated according
     to the scenario tree, we will need to count how many RankCor columns go
     to the diagonal from every period.

     clpddg[0-(level-1)]   number of CORE-columns to go in diagonal from
                           RankCor in period i
     nzpddg[0-(level-1)]   number of nonzeros in above columns
                           (this time sorted by their row-period in CORE)
                          FIXME: is this correct? Do we not need to
                                 distiguish by both row/col period
     nzpddg0               local contribution of current col to nzpddg0
  */

  {
    is_col_diag = Ret->is_col_diag = new int[f_cldiag];

    int *nzpddg0 = new int[nPeriods];
    int *nzpddg  = (int *) calloc(level, sizeof(int));
    int *clpddg  = (int *) calloc(level, sizeof(int));
    bool found;

    // for all columns in border block
    for (i = 0; i < f_cldiag; ++i) {

      // nonzeros in current column
      int nzcl = 0;

      found = false;
      is_col_diag[i] = 0;

      // set cu_pd_cl to correct period block of CORE matrix
      cu_pd_cl = smps.getColPeriod(i);

      // count nonzeros in columns and distribution into row period blocks:
      // - is_col_diag[i]: 1 if no entries in row-blocks assoc with diagonals
      // - nzpddg0[0 - nPeriods-1]: nonzeros for this column in row periods
      // - clpddg[0 - nPeriods-1]: total number of columns that have no entries
      //                           in rows associated with diagonal blocks
      // - nzpddg[0 - nPeriods-1]: total number of entries in above columns
      //                           sorted by which row-periods they occur in

      for (j = 0; j < level; ++j)
	nzpddg0[j] = 0;

      for (j = data.clpnts[i]; j < data.clpnts[i + 1]; ++j) {

	// not objective
	if (data.rwnmbs[j] != objRow) {

	  // get correct row time period for this entry
	  cu_pd_rw = smps.getRowPeriod(data.rwnmbs[j]);

	  // count nonzeros this column will create in final matrix
	  nzpddg0[cu_pd_rw]++;
	  ++nzcl;
	}

	if (data.rwnmbs[j] >= f_rwdiag)
	  found = true;
      }

      // the column can go in the diagonal part
      if (!found) {

	is_col_diag[i] = 1;
	clpddg[cu_pd_cl]++;
	for (j = 0; j < level; ++j)
	  nzpddg[j] += nzpddg0[j];
      }
    }

    // scan through all nodes in first 'level' periods and count
    // col in diagonal and nonzeros in them

    nzdg0 = 0;    // total nonzeros in first diagonal
    cldg0 = 0;    // total columns in first diagonal
    cu_nd_cl = 0;

    node = rootNode;

    do {

      int perNode = node->level();
      if (perNode >= level)
	break;

#ifdef DEBUG
      printf("Node %d period %d,  nz in diag %d, cols in diag %d\n",
	     cu_nd_cl, perNode, nzpddg[perNode], clpddg[perNode]);
#endif
      nzdg0 += nzpddg[perNode];
      cldg0 += clpddg[perNode];
      ++cu_nd_cl;

    } while (node = node->next());

    printf("Cols from RNKCR to DIAG[0]: %d (%d nonzeros).\n", cldg0, nzdg0);

    // clean up
    delete[] nzpddg0;
    free(nzpddg);
    free(clpddg);
  }

  /* nzdg0 = Total nonzeros in final matrix in Border columns that could
             be treated separately (i.e. don't affect periods > level)
     cldg0 = # of columns in final matrix in Border that could be
             treated separately
     cu_nd_cl = total number of nodes in Border */

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

  // count the nonzeros in each period block
  smps.countNzPeriodBlocks();

  // count nonzeros in parts of big matrix
  for (i = 0; i < nPeriods; ++i) {

    rnkc_nz_pd[i] = diag_nz_pd[i] = 0;
    for (j = 0; j < level; ++j)
      rnkc_nz_pd[i] += smps.getNzPeriod(i, j);
    for (j = level; j < nPeriods; ++j)
      diag_nz_pd[i] += smps.getNzPeriod(i, j);
    printf("Nonzeros RNKCR/Diag[%d]: %4d %4d\n",
	   i, rnkc_nz_pd[i], diag_nz_pd[i]);
  }

  node = rootNode;

  do {

    int blkNode = node->block();
    int perNode = node->level();
    rnkc_nz_blk[blkNode] += rnkc_nz_pd[perNode];
    diag_nz_blk[blkNode] += diag_nz_pd[perNode];
    diag_m_blk[blkNode]  += smps.getNRowsPeriod(perNode);
    diag_n_blk[blkNode]  += smps.getNColsPeriod(perNode);
    /*
    printf("%d  %d  %d  ", node->name(), blkNode, perNode);
    printf("* %6d %6d %6d %6d\n",
	   rnkc_nz_blk[blkNode], diag_nz_blk[blkNode],
	   diag_m_blk[blkNode],  diag_n_blk[blkNode]);
    */
  } while (node = node->next());

  rnkc_m_blk[0] = diag_m_blk[0];
  rnkc_n_blk[0] = diag_n_blk[0] - cldg0;
  diag_n_blk[0] = cldg0;

  // save dimensions for postprocess
  Ret->nb_row_rnk  = diag_m_blk[0];
  Ret->nb_col_rnk  = rnkc_n_blk[0];
  Ret->nb_col_diag = diag_n_blk[0];

  // first row/col of given block (in big matrix before reordering)
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
  // allocate sparse the matrices for rankcorrector and diagonal components
  //

  Ret->b = NewDenseVector(ttm, "RHS");
  Ret->c = NewDenseVector(ttn, "Obj");
  Ret->u = NewDenseVector(ttn, "UB");

  // rankcor part
  Algebra **RightColEntries = (Algebra**) malloc((nBlocks + 1) * sizeof(Algebra *));

  // diagonal algebras
  Algebra **DiagEntries = (Algebra**) malloc((nBlocks + 2) * sizeof(Algebra *));

  // bottom algebras (set to zero)
  Algebra **Bottom = (Algebra**) malloc((nBlocks + 1) * sizeof(Algebra *));

  // diagonal elements of Q
  Algebra **QDiag = (Algebra**) malloc((nBlocks + 2) * sizeof(Algebra *));

#ifdef DEBUG
#ifdef WITH_MPI
  if(IS_ROOT_PAR)
#endif
  printf("Dimensions and nonzeros of parts of the deterministic equivalent:\n"
	 " blk |    nonz r/d   |    rows r/d   |    cols r/d   ||"
	 " first row/col\n");
#endif /* DEBUG */

  for (i = 0; i <= nBlocks; ++i) {

#ifdef DEBUG
#ifdef WITH_MPI
  if(IS_ROOT_PAR)
#endif
    printf("%4d | %6d %6d | %6d %6d | %6d %6d || %6d %6d\n", i,
	   rnkc_nz_blk[i], diag_nz_blk[i], rnkc_m_blk[i], diag_m_blk[i],
	   rnkc_n_blk[i],  diag_n_blk[i],  f_rw_blk[i],   f_cl_blk[i]);
#endif /* DEBUG */

    // right-hand columns
    sparse = NewSparseMatrix(rnkc_m_blk[i], rnkc_n_blk[i],
			     rnkc_nz_blk[i], "RankCorPart");
    sparse->cbf = (CallBackFunction) CallBackVoid;
    sparse->nb_el = 0;
    RightColEntries[i] = NewSparseSimpleAlgebra(sparse);

    // diagonal entries
    sparse = NewSparseMatrix(diag_m_blk[i], diag_n_blk[i],
			     diag_nz_blk[i], "DiagPart");
    DiagEntries[i] = NewSparseSimpleAlgebra(sparse);
    sparse->cbf = (CallBackFunction) CallBackVoid;
    sparse->nb_el  = 0;
    sparse->nb_col = 0;

    // bottom part
    sparse = NewSparseMatrix(0, diag_n_blk[i], 0, "BottomPart");
    sparse->cbf = (CallBackFunction) CallBackVoid;
    Bottom[i] = NewSparseSimpleAlgebra(sparse);

    // Q diagonal part
    sparse = NewSparseMatrix(diag_n_blk[i], diag_n_blk[i], 0, "QDiagPart");
    sparse->cbf = (CallBackFunction) CallBackVoid;
    QDiag[i] = NewSparseSimpleAlgebra(sparse);
  }

  sparse = NewSparseMatrix(0, rnkc_n_blk[0], 0, "DiagPart");
  sparse->cbf = (CallBackFunction) CallBackVoid;
  DiagEntries[nBlocks + 1] = NewSparseSimpleAlgebra(sparse);

  sparse = NewSparseMatrix(rnkc_n_blk[0], rnkc_n_blk[0], 0, "QDiagPart");
  sparse->cbf = (CallBackFunction) CallBackVoid;
  QDiag[nBlocks + 1] = NewSparseSimpleAlgebra(sparse);

  //
  // build the deterministic equivalent column by column
  //

  ncol_rc  = 0;     // col added to rankcor so far
  cu_nd_cl = 0;     // node corresponding to current column
  cu_pd_rw = 0;     // period corresponding to current row
  cu_pd_cl = 0;     // period corresponding to current column
  b_cu_blk_cl = 0;  // first column in big matrix of current node
  cu_blk_cl = 0;

  node = rootNode;

  // for all columns in the deterministic equivalent
  for (i = 0; i < ttn; ++i) {

    int perNode = node->level();

    if (i - b_cu_blk_cl >= smps.getNColsPeriod(perNode)) {

      // first col of current node in big matrix
      b_cu_blk_cl = i;

      // current block
      ++cu_nd_cl;
      node = node->next();

      // period of current node
      cu_pd_cl = node->level();

      // current part of big matrix
      cu_blk_cl = node->block();

      /*
      printf("Column %4d - node now %3d  stage now %d  block now %d\n",
	     i, cu_nd_cl, cu_pd_cl, cu_blk_cl);
      */
    }

    // corresponding column in the core matrix
    int coreCol = i - b_cu_blk_cl + smps.getBegPeriodCol(cu_pd_cl);
    assert(coreCol <= smps.getCols());

    // scan through column and set p_pd_rw[pd]:
    // pointers to start of period information in CORE matrix
    // p_pd_rw[pd] is start of info for period 'pd' in current col in CORE

    cu_pd_rw = 0;
    for (k = data.clpnts[coreCol]; k < data.clpnts[coreCol + 1]; ++k) {
      if (k != data.clpnts[coreCol] && (data.rwnmbs[k - 1] > data.rwnmbs[k]))
	printf("WARNING: row numbers are not in ascending order!\n");

      j = data.rwnmbs[k];
      while (j >= smps.getBegPeriodRow(cu_pd_rw)) {
	p_pd_rw[cu_pd_rw] = k;
	++cu_pd_rw;
      }
    }

    while (cu_pd_rw <= nPeriods) {
      p_pd_rw[cu_pd_rw] = data.clpnts[coreCol + 1];
      ++cu_pd_rw;
    }

    // if RankCor Column
    if (cu_blk_cl == 0 && is_col_diag[coreCol] == 0) {

      // initialise col in all RightColEntry matrices
      for (j = 0; j <= nBlocks; ++j) {
	sparse = (SparseSimpleMatrix *) RightColEntries[j]->Matrix;
	sparse->col_beg[ncol_rc] = sparse->nb_el;
	sparse->col_len[ncol_rc] = 0;
      }
      ++ncol_rc;
    }
    // if not RankCor column
    else {
      // initialise col in corresponding diag matrix
      sparse = (SparseSimpleMatrix *) DiagEntries[cu_blk_cl]->Matrix;
      sparse->col_beg[sparse->nb_col] = sparse->nb_el;
      sparse->col_len[sparse->nb_col] = 0;
      sparse->nb_col++;
    }

    // write all parts of this column in the correct RightColEntries,
    // DiagEntries matrices by scanning through the node and its children
    setNodeChildrenRnkc(RightColEntries, DiagEntries,
			p_pd_rw, f_rw_blk, is_col_diag, data,
			node, cu_blk_cl, ncol_rc, coreCol);
  }

  assert(Ret->nb_col_rnk == ncol_rc);

  // write the final pointer for the last+1 column
  for (j = 0; j <= nBlocks; ++j) {
    sparse = (SparseSimpleMatrix *) RightColEntries[j]->Matrix;
    sparse->col_beg[Ret->nb_col_rnk] = sparse->nb_el;

    sparse = (SparseSimpleMatrix *) DiagEntries[j]->Matrix;
    sparse->col_beg[sparse->nb_col] = sparse->nb_el;
  }

  // setup the right-hand side
  setupRhs(smps, Ret);

  // setup objective and bounds
  setupObjective(smps, Ret);

  // apply scenario corrections
  applyScenarios(Ret, DiagEntries, RightColEntries, f_rw_blk, f_cl_blk);

  // reorder objective, bounds and column names
  reorderObjective(Ret, rnkc_n_blk[0]);

  // set up the deterministic equivalent as a DblBordDiagAlgebra
  DblBordDiagSimpleMatrix *MA = NewDblBordDiagSimpleMatrix(nBlocks + 1, Bottom,
							   RightColEntries,
							   DiagEntries,
							   "SMPSAlgA");
  BlockDiagSimpleMatrix *MQ = NewBlockDiagSimpleMatrix(nBlocks + 2, QDiag,
						       "Q_main");

  AlgA = NewDblBordDiagSimpleAlgebra(MA);
  AlgQ = NewBlockDiagSimpleAlgebra(MQ);

  Ret->AlgA = AlgA;
  Ret->AlgQ = AlgQ;

  // clean up
  delete[] p_pd_rw;
  delete[] rnkc_m_blk;
  delete[] rnkc_n_blk;
  delete[] rnkc_nz_blk;
  delete[] rnkc_nz_pd;
  delete[] diag_m_blk;
  delete[] diag_n_blk;
  delete[] diag_nz_blk;
  delete[] diag_nz_pd;
  delete[] f_cl_blk;
  delete[] f_rw_blk;

  return Ret;
}

/** Free the space allocated for the SmpsReturn structure */
void SmpsOops::freeSmpsReturn(SmpsReturn *ret) {

  int i;

  // early return if the structure has not been allocated
  if (!ret)
    return;

  FreeTree(ret->AlgA->Trow);
  FreeTree(ret->AlgA->Tcol);
  FreeTree(ret->AlgQ->Trow);
  FreeTree(ret->AlgQ->Tcol);

  FreeAlgebraAlg(ret->AlgA);
  FreeAlgebraAlg(ret->AlgQ);
  FreeDenseVector(ret->b);
  FreeDenseVector(ret->c);
  FreeDenseVector(ret->u);

  delete[] ret->is_col_diag;

  for (i = 0; i < smps.getTotRows(); ++i)
    delete[] ret->rownames[i];
  delete[] ret->rownames;

  for (i = 0; i < smps.getTotCols(); ++i)
    delete[] ret->colnames[i];
  delete[] ret->colnames;

  delete ret;
}

/**
 * setNodeChildrenRnkc
 *
 * Go through the given node and all its children (recursively): for each
 * node find the correct Block-Matrix (RC[blk], DG[blk]) and append the
 * relevant period part of CORE to it (using the correct row numbers).
 *
 * SETS: - RC[blk], DG[blk], sparse matrix entries
 *
 * @param RC:
 *        Array of rank corrector algebras
 * @param DG:
 *        Array of diagonal algebras
 * @param p_pd_rw:
 *        Pointer to start of this period, this col in CORE
 * @param f_rw_blk[blk]:
 *        First row of BLOCK blk in big matrix
 * @param is_col_diag[i]:
 *        1 if col 'i' should be in Diagonal part;
 *        0 otherwise (Rank Corrector part)
 * @param node:
 *        Node in the deterministic equivalent at which to start
 * @param colBlk:
 *        Block of current column in the deterministic equivalent
 * @param rnkCol:
 *        Index of current column in RankCor part (if in block 0)
 * @param coreCol:
 *        Column in core corresponding to the current column in the
 *        deterministic equivalent
 */
void SmpsOops::setNodeChildrenRnkc(Algebra **RC, Algebra **DG,
				   int *p_pd_rw, int *f_rw_blk, int *is_col_diag,
				   const SparseData &data,
				   const Node *node, const int colBlk,
				   const int rnkCol, const int coreCol) {

  const int per = node->level();
  const int blk = node->block();
  int index;
  SparseSimpleMatrix *sparse;

  // rank corrector
  if (colBlk == 0 && is_col_diag[coreCol] == 0) {

    sparse = (SparseSimpleMatrix *) RC[blk]->Matrix;
    index  = rnkCol;
  }

  // diagonal
  else {
    if (blk != colBlk && p_pd_rw[per] != p_pd_rw[per + 1]) {
      printf("Entry in non diagonal part of diagonal\n");
      exit(1);
    }

    sparse = (SparseSimpleMatrix *) DG[blk]->Matrix;
    index  = sparse->nb_col;
  }

  // copy the information for this node into the deterministic equivalent
  for (int k = p_pd_rw[per]; k < p_pd_rw[per + 1]; ++k) {

    sparse->element[sparse->nb_el] = data.acoeff[k];
    sparse->row_nbs[sparse->nb_el] = data.rwnmbs[k]
      - smps.getBegPeriodRow(per) + node->firstRow()
      - f_rw_blk[blk];

    /*
    printf(" %2d  - %2d :> %2d, %2d, %2d, %2d  ", per, ordNode,
	   data.rwnmbs[k], smps.getBegPeriodRow(per),
	   node->firstRow(), f_rw_blk[blk]);
    printf(":: %lf  %d\n", sparse->element[sparse->nb_el],
	   sparse->row_nbs[sparse->nb_el]);
    */

    sparse->nb_el++;
    sparse->col_len[index - 1]++;

    assert(sparse->nb_el <= sparse->max_nb_el);
  }

  for (int i = 0; i < node->nChildren(); ++i) {
    setNodeChildrenRnkc(RC, DG, p_pd_rw, f_rw_blk, is_col_diag,
			data, node->getChild(i), colBlk, rnkCol, coreCol);
  }
}

/** Set up the right-hand side */
void setupRhs(const Smps &smps, SmpsReturn *Ret) {

  int firstRowNode, begRowPeriod, period;
  char scname[8], *p;
  DenseVector *rhs = Ret->b;
  const Node *node = Ret->rootNode;

  // leave immediately if there is no root node
  if (!node)
    return;

  char **rownames = Ret->rownames = new char*[smps.getTotRows()];

  // for all nodes in the tree in order
  do {

    period = node->level();
    firstRowNode = node->firstRow();
    begRowPeriod = smps.getBegPeriodRow(period);

    sprintf(scname, "_S%03d", node->name());

    // copy the information for this node
    for (int i = 0; i < node->nRows(); ++i) {

      rhs->elts[firstRowNode + i] = smps.getRhs(begRowPeriod + i);

      // build a name for this row
      rownames[firstRowNode + i] = new char[15];
      strncpy(rownames[firstRowNode + i], smps.getBegRowName(begRowPeriod + i), 8);
      p = rownames[firstRowNode + i] + 8;
      while(*(p - 1) == ' ') --p;
      strcpy(p, scname);
      /*
      printf("%s\n", rownames[firstRowNode + i]);
      */
    }

  } while (node = node->next());
}

/** Set up the objective and the bounds */
void setupObjective(const Smps &smps, SmpsReturn *Ret) {

  int row, firstColNode, begColPeriod, period;
  char buffer[50], scname[8], *p;
  DenseVector *obj = Ret->c, *upb = Ret->u;
  const Node *node = Ret->rootNode;

  // leave immediately if there is no root node
  if (!node)
    return;

  char **colnames = Ret->colnames = new char*[smps.getTotCols()];

  // copy the objective row from the core matrix
  double *coreObj = smps.getObjRow();

  // for all nodes in the tree
  do {

    period  = node->level();
    firstColNode = node->firstCol();
    begColPeriod = smps.getBegPeriodCol(period);
    double probNode = node->probNode();

    sprintf(scname, "_S%03d", node->name());

    // copy the information for this node
    for (int i = 0; i < node->nCols(); ++i) {

      // copy the objective coefficients weighted by probability of the node
      obj->elts[firstColNode + i] = probNode * coreObj[begColPeriod + i];

      // copy the upper bounds
      upb->elts[firstColNode + i] = smps.getUpperBound(begColPeriod + i);

      // and build a name for this col
      colnames[firstColNode + i] = new char[20];
      strncpy(buffer, smps.getBegColName(begColPeriod + i), 8);
      p = buffer + 8;
      if (strncmp(buffer, "SK", 2) == 0) {
	row = atoi(buffer + 2);
	buffer[2] = '_';
	strncpy(buffer + 3, smps.getBegRowName(row), 8);
	p = buffer + 11;
      }
      while(*(p - 1) == ' ') --p;
      *p = 0;
      strcpy(colnames[firstColNode + i], buffer);
      strcat(colnames[firstColNode + i], scname);
      /*
      printf("%s\n",colnames[firstColNode+i]);
      */
    }
  } while (node = node->next());

  // clean up
  delete[] coreObj;
}

/**
 *  Apply the scenarios changes.
 *
  IN:
   scenario change data:
     sc_first[sc]   pointer to first entry for scenario sc
     sc_len[sc]     number of entries for scenario sc
     row_nbs[i]     affected row
     col_nbs[i]     affected col
     entry[i]       new entry
     NOTE: sc_first, row_nbs, col_nbs are in FORTRAN numbnering!
   obj,rhs,matrix:
     Ret->b           current RHS
     Ret->c           current objective
                   these two are arrays of double in a depth first node
                   ordering. Starting rows and columns corresponding to each
                   node are given by f_cl_nd, f_rw_nd in SmpsTree.
                   These vectors are reordered later.
     DiagEntries      Diagonal Matrix blocks
     RightColEntries  Right Column Matrix blocks
                   these are arrays of Algebras. They are already in final
		   ordering. The ordering is as follows
                   Variables:
                     RnkCor: proper first stage decision variable
                     D[0]  : first stage variables that are not refered to by
                             later stages
                   Constraints:
                     D[1] -D[n]: scenario blocks
                     D[0]      : internal first stage constraints

  Need mapping from DetEq row/col to Block number in
  DiagEntries/RightColEntries and row/col number within these.
  Within each matrix blocks nodes are ordered in a depth first ordering

  block[i]:          number of block that node belongs to
  scenario[i]        scenario that node belongs to
  period[i]          period that node belongs to

  Need: for each node, start of row/col within its DetEquivMatrix blck
    -   for diagonal blocks can take row/col from f_cl_nd, f_rw_nd and
        substract starting row/col for seed node of this block
    -   for nondiag node need to do breadth first count
        the correct column within the RncD0 part is obtained
	by looping through all columns within RncD0 (breadth first over nodes)
	and counting how many end up in Rnc and D0
  OUT: changed RHSm objective, system matrix

  Every scenario is a path from the root node to one of the leaf
  nodes.  the scenario change data has for all scenarios a list of
  changes that need to be applied. These changes can apply to any
  period of the tree (not just the leaf nodes).

  FIXME: we should do this by looping over scenarios first

  We apply the scenario change data by the following algorithm

  LOOP over all leaf nodes lf

     LOOP walk up from leaf node to root node: scNode

       work out which scenario scNode corresponds to: scen

       scan through all the scenario changes:
         Scenario changes are stored by entries in rhs,obj,matrix of the core

         get row/col period of this change

         apply only changes that correspond to this period

         if (objective change)
           work out column of change and apply
         if (rhs change)
           work out row of change and apply
	 if (matrix change)
	   work out row and column of change and apply

     END

  END
*/
int SmpsOops::applyScenarios(SmpsReturn *Ret,
			     Algebra **DiagEntries, Algebra **RightColEntries,
			     int *f_rw_blk, int *f_cl_blk) {

  int period, firstEntry, lastEntry;
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

      // scenario the current node belongs to
      int scen = scNode->scenario();

      // index of the first change
      firstEntry = smps.getFirstEntryScen(scen) - 1;
      lastEntry  = firstEntry + smps.getLengthScen(scen);

#ifdef DEBUG_SCEN
      printf("Node %d: scenario %d (entries from %d to %d)\n",
	     scNode->name(), scen + 1, firstEntry, lastEntry - 1);
#endif

      // for all changes affecting this scenario
      for (int corr = firstEntry; corr < lastEntry; ++corr) {

	// row and column of core affected by the change
	pdr = smps.getRowPeriod(entryRow[corr] - 1);
	pdc = smps.getColPeriod(entryCol[corr] - 1);

	// if the change affects the objective
	if (pdr < 0 && pdc == period) {

	  col = scNode->firstCol() + entryCol[corr] - 1
	    - smps.getBegPeriodCol(pdc);
	  assert(col <= smps.getTotCols());

	  Ret->c->elts[col] = entryVal[corr] * scNode->probNode();

#ifdef DEBUG_SCEN
	  printf("   Obj entry: core col %d, det.eq. col %d, (%f)\n",
		 entryCol[corr], col, entryVal[corr]);
#endif
	}

	// if the change affects the rhs
	else if (pdc < 0 && pdr == period) {

	  row = scNode->firstRow() + entryRow[corr] - 1
	    - smps.getBegPeriodRow(pdr);
	  assert(row <= smps.getTotRows());

	  Ret->b->elts[row] = entryVal[corr];

#ifdef DEBUG_SCEN
	  printf("   Rhs entry: core row %d, det.eq. row %d, (%f)\n",
		 entryRow[corr], row, entryVal[corr]);
#endif
	}

	// if the change affects the matrix
	else if (pdc >= 0 && pdr >= 0 && pdr == period) {

	  // indices in the deterministic equivalent
	  row = scNode->firstRow() + entryRow[corr] - 1
	    - smps.getBegPeriodRow(pdr);
	  col = scNode->firstCol() + entryCol[corr] - 1
	    - smps.getBegPeriodCol(pdc);

	  // adjust the column index if the change is in a column that belongs
	  // to the previous period
	  if (pdr != pdc)
	    col += scNode->parent()->firstCol() - scNode->firstCol();

	  assert(row <= smps.getTotRows());
	  assert(col <= smps.getTotCols());

#ifdef DEBUG_SCEN
	  printf("   Mtx entry: core row %3d col %3d, det.eq. row %3d col %3d",
		 entryRow[corr], entryCol[corr], row, col);
#endif

	  // work out the correct row and column block of the correction
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
	      sparse = (SparseSimpleMatrix *) DiagEntries[0]->Matrix;
	      iblk = iblkd0;
	    } else {
	      sparse = (SparseSimpleMatrix *) RightColEntries[rowBlock]->Matrix;
	      iblk = iblkrnc;
	    }
	  }

	  // if the change is in the diagonal part
	  else {

	    // ensure we are in a diagonal block
	    assert(colBlock == rowBlock);

	    sparse = (SparseSimpleMatrix *) DiagEntries[colBlock]->Matrix;
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
	    printf("Entry not found\n");
	    exit(1);
	  }
	}
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
void SmpsOops::reorderObjective(SmpsReturn *Ret, const int rnkn) {

  int col, coreCol, firstColDiag, firstColNode;
  int nb_el = 0;
  int ttn = smps.getTotCols();

  const Node *node = Ret->rootNode;
  if (!node)
    return;

  DenseVector *obj = Ret->c, *upb = Ret->u;
  int *is_col_diag = Ret->is_col_diag;
  char **colnames  = Ret->colnames;

  double *objCopy = (double *) malloc(ttn * sizeof(double));
  double *upbCopy = (double *) malloc(ttn * sizeof(double));
  char  **clnCopy = (char **)  malloc(ttn * sizeof(char *));

  // find the first node that will go in a diagonal block
  while (node->level() < level)
    node = node->next();

  // set the first column in diagonal block
  firstColDiag = node->firstCol();

  // copy objective and bounds into temporary arrays
  for (int i = 0; i < ttn; ++i) {
    objCopy[i] = obj->elts[i];
    upbCopy[i] = upb->elts[i];
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
      upb->elts[nb_el] = upbCopy[col];
      colnames[nb_el]  = clnCopy[col];
      ++nb_el;
    }
  }

  assert(nb_el == firstColDiag - rnkn);

  // copy the other diagonal entries  (Diag-1 ... Diag-n)
  for (col = firstColDiag; col < ttn; ++col) {
    obj->elts[nb_el] = objCopy[col];
    upb->elts[nb_el] = upbCopy[col];
    colnames[nb_el]  = clnCopy[col];
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
      upb->elts[nb_el] = upbCopy[col];
      colnames[nb_el]  = clnCopy[col];
      ++nb_el;
    }
  }

  assert(nb_el == ttn);

  // clean up
  free(objCopy);
  free(upbCopy);
  free(clnCopy);
}

/**
 *  Copy a Vector into a depth-first ordered DenseVector.
 *
 *  Copies a Vector as used by OOPS into a DenseVector corresponding to a
 *  depth-first ordering of the scenario tree.
 */
void SmpsOops::SmpsVectorToDense(Vector *x, DenseVector *dx,
				 SmpsReturn *Ret, const int rowcol) {

  SetExactVector(x);
  CopyToDenseVector(x, dx);

  // reorder according to the SMPS depth-first order
  if (rowcol == ORDER_COL)
    backOrderColVector(dx->elts, Ret);
  else
    backOrderRowVector(dx->elts, Ret);
}

/**
 *  Copy a depth-first ordered DenseVector into a Vector.
 *
 *  Copies a DenseVector corresponding to a depth-first ordering of the
 *  scenario tree into a Vector as used by OOPS.
 *
 *  Reorders the elements according to the reordering used by OOPS:
 *  first period entries are placed at the end rather than at the beginning,
 *  those columns that are not linking periods are placed in separate
 *  diagonal block, rather than in the RankCor block.
 */
void SmpsOops::SmpsDenseToVector(DenseVector *dx, Vector *x,
				 SmpsReturn *Ret, const int rowcol) {

  // should attempt to leave the original element order intact

  // go for memory saving option
  if (rowcol == ORDER_COL)
    forwOrderColVector(dx->elts, Ret);
  else
    forwOrderRowVector(dx->elts, Ret);

  CopyDenseToVector(dx, x);

  // and reverse the order, to leave the dense vector intact
  if (rowcol == ORDER_COL)
    backOrderColVector(dx->elts, Ret);
  else
    backOrderRowVector(dx->elts, Ret);
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
 *  - nb_col_rnk:  number of columns in actual rankcor (Rnk) part
 *  - nb_col_diag: number of columns in diag rankcor (D0) part
 */
void SmpsOops::backOrderColVector(double *x, const SmpsReturn *Ret) {

  const Node *node = Ret->rootNode;

  // total number of columns in RankCor (D0|Rnk)
  const int ncol_ttrc = Ret->nb_col_rnk + Ret->nb_col_diag;
  const int ttn = smps.getTotCols();
  double *dtmp  = new double[ttn];

  for (int i = 0; i < ttn; ++i)
    dtmp[i] = x[i];

  // copy entries into the combined rankcor slot

  // set the start of the D0 and Rnk blocks
  int nx_col_d0 = 0;                      // next column to take from D0
  int nx_col_rc = ttn - Ret->nb_col_rnk;  // next column to take from Rnk
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
    if (Ret->is_col_diag[coreCol])
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
 *  - nb_col_rnk:  number of columns in actual rankcor (Rnk) part
 *  - nb_col_diag: number of columns in diag rankcor (D0) part
 */
void SmpsOops::forwOrderColVector(double *x, const SmpsReturn *Ret) {

  const Node *node = Ret->rootNode;

  // total number of columns in RankCor (D0|Rnk)
  const int ncol_ttrc = Ret->nb_col_rnk + Ret->nb_col_diag;
  const int ttn = smps.getTotCols();
  double *dtmp  = new double[ttn];

  for (int i = 0; i < ttn; ++i)
    dtmp[i] = x[i];

  // copy entries from the combined rankcor slot into OOPS RankCor and
  // first diagonal

  // set the start of the D0 and Rnk blocks
  int nx_col_d0 = 0;                      // next column in D0
  int nx_col_rc = ttn - Ret->nb_col_rnk;  // next column in Rnk
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
    if (Ret->is_col_diag[coreCol])
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

/** Reorder the rows */
void SmpsOops::backOrderRowVector(double *x, const SmpsReturn *Ret) {

  int i;
  const int ttm = smps.getTotRows(), nRowsRnkc = Ret->nb_row_rnk;
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
 *  Order the rows from the SMPS depth-first order into the order that is used
 *  internally by OOPS (the only difference is that OOPS has the rows
 *  correspoding to the first 'level' periods at the end, whereas the SMPS
 *  depth-first order has them at the beginning).
 *
 *  @param x:
 *         The vector that should be reordered
 *  @param Ret:
 *         Information about the problem with respect to which the
 *         reordering should be done
 */
void SmpsOops::forwOrderRowVector(double *x, const SmpsReturn *Ret) {

  int i;
  const int ttm = smps.getTotRows(), nRowsRnkc = Ret->nb_row_rnk;
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
