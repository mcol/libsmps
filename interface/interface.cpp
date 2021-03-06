/*
 *  interface.cpp
 *
 *  Generate the deterministic equivalent from the SMPS data structures.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "interface.h"

static ProbData* setupMatrix(Smps &smps, const Node *node);
static int setupRhs(ProbData *PROB, const Smps &smps, const Node *node);
static int applyScenarios(ProbData *PROB, const Smps &smps, const Node *node);
static int copyLinkingBlocks(Smps &smps, SparseData &data, const Node *node,
			     double *acoeff, int *rwnmbs,
			     int &cIndex, int &dIndex, int &nnzCol);


/** Generate the deterministic equivalent problem */
ProbData* setupProblem(Smps &smps, const Node *node) {

  ProbData *PROB;

  printf(" --------------- setupProblem --------------\n");

  // setup the matrix
  PROB = setupMatrix(smps, node);
  if (!PROB)
    return NULL;

  // setup the right-hand side
  setupRhs(PROB, smps, node);

  // apply the scenario corrections
  applyScenarios(PROB, smps, node);

  return PROB;
}

/** Setup the constraint matrix */
ProbData *setupMatrix(Smps &smps, const Node *node) {

  if (!node)
    node = smps.getRootNode();

  // leave immediately if there is no root node
  if (!node)
    return NULL;

  int period, col = 0;
  int curCol, nnzCol;
  int f_cl_nd, f_cl_pd;
  int *rwnmbs, *clpnts;
  double *acoeff, *obj, *blo, *bup;

  int cIndex = 0;   // index in the core matrix
  int dIndex = 0;   // index in the deterministic equivalent

  int iob = smps.getObjRowIndex();
  int ttm = smps.getTotRows();
  int ttn = smps.getTotCols();

  // these assertions may be triggered if the time file leaves
  // some rows or columns unassigned to any period
  assert(ttm >= smps.getRows() - 1);
  assert(ttn >= smps.getCols());

  // count the number of nonzero elements
  int ttnz = smps.countNonzeros(smps.getSmpsTree());

  // dimensions of the deterministic equivalent
  printf("The deterministic equivalent matrix is %dx%d, %d nonzeros.\n",
	 ttm, ttn, ttnz);

  // allocate space for the arrays
  obj    = new double[ttn];
  blo    = new double[ttn];
  bup    = new double[ttn];
  acoeff = new double[ttnz];
  rwnmbs = new int[ttnz];
  clpnts = new int[ttn + 1];
  if (obj    == NULL || blo    == NULL || bup    == NULL ||
      acoeff == NULL || rwnmbs == NULL || clpnts == NULL) {
    printf("Memory allocation failed.\n");
    delete[] obj;
    delete[] blo;
    delete[] bup;
    delete[] acoeff;
    delete[] rwnmbs;
    delete[] clpnts;
    return NULL;
  }
  memset(blo, 0, ttn * sizeof(double));

  SparseData data = smps.getSparseData();

  // for all nodes
  do {

    // find which period this node belongs to
    period = node->level();

    // offset in row numbers between core and deterministic equivalent
    int offset = node->firstRow() - smps.getBegPeriodRow(period);

    // first column of the period for the deterministic equivalent and core
    f_cl_nd = node->firstCol();
    f_cl_pd = smps.getBegPeriodCol(period);

    // initialize objective and bounds
    for (int i = 0; i < node->nCols(); ++i) {

      obj[f_cl_nd + i] = 0.0;
      blo[f_cl_nd + i] = smps.getLowerBound(f_cl_pd + i);
      bup[f_cl_nd + i] = smps.getUpperBound(f_cl_pd + i);
    }

    // index in core corresponding to the beginning of this period
    cIndex = data.clpnts[f_cl_pd];

    // for all columns in this period
    while (cIndex < data.clpnts[smps.getBegPeriodCol(period + 1)]) {

      assert(col <= ttn);
      assert(dIndex < ttnz);

      // set the column pointer
      clpnts[col] = dIndex;

      // index of the current column in core
      curCol = data.clnmbs[cIndex];

      // number of nonzeros in this column
      nnzCol = data.clpnts[curCol + 1] - data.clpnts[curCol];

      // copy the rows of the period
      while (nnzCol > 0) {

	int row = data.rwnmbs[cIndex];

	// stop early before starting a new period
	if (row >= smps.getBegPeriodRow(period + 1))
      	  break;

	// objective row
	if (row == iob) {

	  obj[col] = data.acoeff[cIndex] * node->probNode();
	}

	// this element is above the diagonal
	else if (row < smps.getBegPeriodRow(period)) {

	  int pdr = smps.getRowPeriod(row);
	  int pdc = smps.getColPeriod(curCol);
	  assert(pdc > pdr);

	  // node corresponds to the column period
	  // find the node corresponding to the row period
	  const Node *rowNode = node;
	  while (rowNode->level() > pdr)
	    rowNode = rowNode->parent();
	  rwnmbs[dIndex] = row
	    + rowNode->firstRow() - smps.getBegPeriodRow(pdr);

	  // scale the element by the conditional probability
	  acoeff[dIndex] = data.acoeff[cIndex] 
	    * node->probNode() / rowNode->probNode();

#ifdef DEBUG_MATRIX
	  printf("Got it! %3g, %d (per %d) node %d\t",
		 data.acoeff[cIndex], row, period, node->name());
	  printf("pdr: %d, pdc: %d\n", pdr, pdc);
	  printf("%2d (%2d)| %10f  place at (%2d, %2d)\n",
		 dIndex, cIndex, acoeff[dIndex], rwnmbs[dIndex], col);
#endif
	  assert(rwnmbs[dIndex] >= 0);
	  assert(rwnmbs[dIndex] < ttm);

	  ++dIndex;
	}

	else {
	  acoeff[dIndex] = data.acoeff[cIndex];
	  rwnmbs[dIndex] = row + offset;

	  assert(rwnmbs[dIndex] >= 0);
	  assert(rwnmbs[dIndex] < ttm);

#ifdef DEBUG_MATRIX
	  printf("%2d (%2d)| %10f  %d\n",
		 dIndex, cIndex, acoeff[dIndex], rwnmbs[dIndex]);
#endif
	  ++dIndex;
	}

	++cIndex;
	--nnzCol;
      }

      // copy the linking blocks in the current period
      copyLinkingBlocks(smps, data, node,
			acoeff, rwnmbs, cIndex, dIndex, nnzCol);

      // go to the next column
      ++col;

    } // end while

  } while ((node = node->next()));

  assert(col == ttn);

  // set the last column pointer
  clpnts[col] = dIndex;

  // allocate the ProbData structure
  ProbData *PROB = new ProbData;

  // fill the ProbData structure
  PROB->ttm    = ttm;
  PROB->ttn    = ttn;
  PROB->ttnz   = clpnts[ttn];
  PROB->irobj  = iob;
  PROB->acoeff = acoeff;
  PROB->rwnmbs = rwnmbs;
  PROB->clpnts = clpnts;
  PROB->obj    = obj;
  PROB->blo    = blo;
  PROB->bup    = bup;
  PROB->rwnames = smps.getRowNames(smps.getSmpsTree());
  PROB->clnames = smps.getColNames(smps.getSmpsTree());

  return PROB;
}

/** Copy the elements in the linking blocks in the current column */
int copyLinkingBlocks(Smps &smps, SparseData &data, const Node *node,
		      double *acoeff, int *rwnmbs,
		      int &cIndex, int &dIndex, int &nnzCol) {

  const int ttm = smps.getTotRows();

  // store the current row index of core
  const int sIndex = cIndex;

  // store the number of remaining nonzeros in the column
  const int snzCol = nnzCol;

  // the period of the node
  const int period = node->level();

  // copy the linking blocks from the current period
  for (int block = 0; block < node->nChildren(); ++block) {

    const Node *child = node->getChild(block);

    // determine the offset in row numbers for the children nodes
    int offsetChild = child->firstRow() - smps.getBegPeriodRow(period + 1);

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

      // this element belongs to the next period
      if (data.rwnmbs[cIndex] < smps.getBegPeriodRow(period + 2)) {
	acoeff[dIndex] = data.acoeff[cIndex];
	rwnmbs[dIndex] = data.rwnmbs[cIndex] + offsetChild;
#ifdef DEBUG_MATRIX
	printf("%2d (%2d)| %10f  %d from node %d\n",
	       dIndex, cIndex, acoeff[dIndex], rwnmbs[dIndex], node->name());
#endif

	assert(rwnmbs[dIndex] >= 0);
	assert(rwnmbs[dIndex] < ttm);

	++cIndex;
	++dIndex;
	--nnzCol;
      }

      // this element belongs to a period after the next
      else
	copyLinkingBlocks(smps, data, child,
			  acoeff, rwnmbs, cIndex, dIndex, nnzCol);
    }
  }

  return 0;
}

/** Setup the right-hand side */
int setupRhs(ProbData *PROB, const Smps &smps, const Node *node) {

  if (!node)
    node = smps.getRootNode();

  // leave immediately if there is no root node
  if (!node)
    return 1;

  // allocate space for right-hand side and row types
  double *rhs = new double[PROB->ttm];
  int    *rws = new int[PROB->ttm];

  // for all nodes
  do {

    // find which period the node belongs to
    int period = node->level();

    // index of the first row
    int firstRow = node->firstRow();

    // for all rows in this period
    for (int i = 0; i < node->nRows(); ++i) {

      int index = firstRow + i;
      assert (index < PROB->ttm);

      rhs[index] = smps.getRhs(smps.getBegPeriodRow(period) + i);
      rws[index] = smps.getRowType(smps.getBegPeriodRow(period) + i);
    }

  }  while ((node = node->next()));

  PROB->rhs = rhs;
  PROB->rws = rws;

  return 0;
}

/** Apply the scenario corrections */
int applyScenarios(ProbData *PROB, const Smps &smps, const Node *node) {

  int period, firstEntry, lastEntry;
  int row, col, pdr, pdc, index;
  const int *entryRow = smps.getEntryRow();
  const int *entryCol = smps.getEntryCol();
  const double *entryVal = smps.getEntryVal();
  const Node *scNode;

  if (!node)
    node = smps.getRootNode();

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

	  PROB->obj[col] = entryVal[corr] * scNode->probNode();

#ifdef DEBUG_SCEN
	  printf("   Obj entry: core col %3d, det.eq. col %3d, (%f)\n",
		 entryCol[corr], col, entryVal[corr]);
#endif
	}

	// if the change affects the rhs
	else if (pdc < 0 && pdr == period) {

	  row = scNode->firstRow() + entryRow[corr] - 1
	    - smps.getBegPeriodRow(pdr);
	  assert(row <= smps.getTotRows());

	  PROB->rhs[row] = entryVal[corr];

#ifdef DEBUG_SCEN
	  printf("   Rhs entry: core row %3d, det.eq. row %3d, (%f)\n",
		 entryRow[corr], row, entryVal[corr]);
#endif
	}

	// the change is in an above-diagonal element
	else if (pdr >= 0 && pdc > pdr && pdc == period) {

	  // node corresponds to the column period
	  // find the node corresponding to the row period
	  const Node *rowNode = scNode;
	  while (rowNode->level() > pdr)
	    rowNode = rowNode->parent();

	  // indices in the deterministic equivalent
	  row = rowNode->firstRow() + entryRow[corr] - 1
	    - smps.getBegPeriodRow(pdr);
	  col = scNode->firstCol() + entryCol[corr] - 1
	    - smps.getBegPeriodCol(pdc);

	  assert(row <= smps.getTotRows());
	  assert(col <= smps.getTotCols());

	  index = PROB->clpnts[col];
	  while (PROB->rwnmbs[index] != row)
	    index++;

#ifdef DEBUG_SCEN
	  printf("   Abd entry: core row %3d col %3d, det.eq. row %3d col %3d,"
		 " (%g -> %g)\n",
		 entryRow[corr], entryCol[corr], row, col,
		 PROB->acoeff[index], entryVal[corr]);
#endif
	  PROB->acoeff[index] = entryVal[corr]
	    * scNode->probNode() / rowNode->probNode();
	}

	// if the change affects the matrix
	else if (pdc >= 0 && pdr >= pdc && pdr == period) {

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

	  assert(row <= smps.getTotRows());
	  assert(col <= smps.getTotCols());

	  index = PROB->clpnts[col];
	  while (PROB->rwnmbs[index] != row)
	    index++;

#ifdef DEBUG_SCEN
	  printf("   Mtx entry: core row %3d col %3d, det.eq. row %3d col %3d,"
		 " (%g -> %g)\n",
		 entryRow[corr], entryCol[corr], row, col,
		 PROB->acoeff[index], entryVal[corr]);
#endif
	  PROB->acoeff[index] = entryVal[corr];
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
      --period;
    }
  } while ((node = node->next()));

  return 0;
}

/** Free the space allocated for the ProbData structure */
int freeProbData(ProbData *PROB) {

  delete[] PROB->acoeff;
  delete[] PROB->rwnmbs;
  delete[] PROB->clpnts;
  delete[] PROB->obj;
  delete[] PROB->blo;
  delete[] PROB->bup;
  delete[] PROB->rhs;
  delete[] PROB->rws;

  if (PROB->rwnames) {
    for (int i = 0; i < PROB->ttm; ++i)
      delete[] PROB->rwnames[i];
    delete[] PROB->rwnames;
  }

  if (PROB->clnames) {
    for (int i = 0; i < PROB->ttn; ++i)
      delete[] PROB->clnames[i];
    delete[] PROB->clnames;
  }

  delete PROB;

  return 0;
}

/** Print the solution */
void printSolution(const Node *root,
		   double *primal, double *dual,
		   double *slacks, double *rcosts) {

  int i, end;
  const Node *node = root;

  // return immediately if there is no root node
  if (!node)
    return;

  // open the output file
  FILE *out = fopen("smps.sol", "w");

  // print the solution in breadth-first order
  queue<const Node*> qNodes;

  qNodes.push(node);

  while (!qNodes.empty()) {

    node = qNodes.front();
    qNodes.pop();

    fprintf(out, "\t---   Node %2d (%dx%d)  ---\n", node->name(),
            node->nRows(), node->nCols());

    end = node->firstRow() + node->nRows();
    for (i = node->firstRow(); i < end; ++i) {
      fprintf(out, "Row %d:  Slack = %10f  Dual = %10f\n", i,
              slacks[i], dual[i]);
    }

    end = node->firstCol() + node->nCols();
    for (i = node->firstCol(); i < end; ++i) {
      fprintf(out, "Column %d:  Value = %10f  Reduced cost = %10f\n", i,
              primal[i], rcosts[i]);
    }

    // add the children to the queue of nodes to print
    for (i = 0; i < node->nChildren(); ++i)
      qNodes.push(node->getChild(i));
  }

  // close the output file
  fclose(out);
}

#ifdef OBSOLETE
/**
 *  Write the deterministic equivalent problem in MPS format.
 *
 *  MPS format for QP is as follows:
 *  000 0    1    1    2    2    3    3    4    4    5    5    66
 *  123 5    0    5    0    5    0    5    0    5    0    5    01
 *  -------------------------------------------------------------
 *  NAME          [name of problem]
 *  ROWS
 *   E  [name]          (equality "=" c/s)
 *   G  [name]          (greater than or equal ">=" c/s)
 *   L  [name]          (less than or equal "<=" c/s)
 *   N  [name]          (first N row found is objective)
 *  COLUMNS
 *      COL1      ROW1      1.0       ROW2      -1.0
 *      COL1      ROW3      4.0
 *      COL2      ROW1      2.3       ROW3      0.4
 *       ....
 *   (gives all nonzero entries of A, column by column)
 *   (giving a second entry in same line is optional)
 *  RHS
 *      rhs       ROW1      3.0       ROW2      4,9
 *      rhs       ROW3      8.9
 *   (all other rhs entries assumed to be 0)
 *  RANGES
 *   (optional allows upper and lower bounds on constraints to be defined)
 *  BOUNDS
 *   LO bnd       COL1      0.4
 *   UB bnd       COL1      0.6
 *   (LO defines lower bound, UB defines upper bound)
 *  ENDATA
 */
int writeMps(const char *filename, ProbData *PROB) {

  /* format of the mps strings */
  const char cl_row_format[] = " %c  R%-7d\n";
  const char cl_obj_format[] = "    C%-7d  OBJ       %-10.7G\n";
  const char cl_col_format[] = "    C%-7d  R%-7d  %-10.7G\n";
  const char cl_rhs_format[] = "    RHS       R%-7d  %-10.7G\n";
  const char cl_lob_format[] = " LO BOUND     C%-6d    %10.7G\n";
  const char cl_upb_format[] = " UP BOUND     C%-6d    %10.7G\n";

  const char* row_format = cl_row_format;
  const char* obj_format = cl_obj_format;
  const char* col_format = cl_col_format;
  const char* rhs_format = cl_rhs_format;
  const char* lob_format = cl_lob_format;
  const char* upb_format = cl_upb_format;

  /* convert PROB->rws into the corresponding letter */
  const char rowtype[] = "?EGLN";

  /* small value */
  const double epsilon = 1.E-30;

  int row, col;

  /* open the output file for writing */
  FILE *out = fopen(filename, "w");
  if (!out)
    return 1;

  /* NAME section */
  fprintf(out, "NAME SmpsGen\n");

  /* ROWS section */
  fprintf(out, "ROWS\n");
  fprintf(out, " N  OBJ\n");
  for (row = 0; row < PROB->ttm; ++row)
    fprintf(out, row_format, rowtype[PROB->rws[row]], row);

  /* COLUMNS section */
  fprintf(out, "COLUMNS\n");
  for (col = 0; col < PROB->ttn; ++col) {

    int j;

    /* write the objective row */
    if (PROB->obj[col] > epsilon || PROB->obj[col] < -epsilon)
      fprintf(out, obj_format, col, PROB->obj[col]);

    /* write the other rows */
    for (j = PROB->clpnts[col]; j < PROB->clpnts[col + 1]; ++j) {

      row = PROB->rwnmbs[j];
      fprintf(out, col_format, col, row, PROB->acoeff[j]);
    }
  }

  /* RHS section */
  fprintf(out, "RHS\n");
  for (row = 0; row < PROB->ttm; ++row)
    if (PROB->rhs[row] > 0.)
      fprintf(out, rhs_format, row, PROB->rhs[row]);

  /* BOUNDS section */
  fprintf(out, "BOUNDS\n");
  for (col = 0; col < PROB->ttn; ++col) {
    if (PROB->blo[col] > 0.)
      fprintf(out, lob_format, col, PROB->blo[col]);
    if (PROB->bup[col] < 1e31)
      fprintf(out, upb_format, col, PROB->bup[col]);
  }

  /* end of the file */
  fprintf(out, "ENDATA\n");

  /* close the output file */
  fclose(out);
  printf("Deterministic equivalent matrix written in file '%s'.\n", filename);

  return 0;
}
#endif /* OBSOLETE */
