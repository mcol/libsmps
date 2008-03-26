/*
 *  setupProblem.c
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
#include <assert.h>
#include "interface.h"

static ProbData* setupMatrix(Smps &smps);
static int setupRhs(ProbData *PROB, const Smps &smps);
static int applyScenarios(ProbData *PROB, const Smps &smps);


/** Generate the deterministic equivalent problem */
ProbData* setupProblem(Smps &smps) {

  ProbData *PROB;

  printf(" --------------- setupProblem --------------\n");

  // setup the matrix
  PROB = setupMatrix(smps);
  if (!PROB)
    return NULL;

  // setup the right-hand side
  setupRhs(PROB, smps);

  // apply the scenario corrections
  applyScenarios(PROB, smps);

  return PROB;
}

/** Setup the constraint matrix */
ProbData *setupMatrix(Smps &smps) {

  const Node *node = smps.getRootNode();

  // leave immediately if there is no root node
  if (!node)
    return NULL;

  int period, col = 0;
  int curCol, nnzCol, snzCol, sIndex;
  int f_cl_nd, f_cl_pd;
  int *rwnmbs, *clpnts;
  double *acoeff, *obj, *blo, *bup;

  int cIndex = 0;   // index in the core matrix
  int dIndex = 0;   // index in the deterministic equivalent

  int iob = smps.getObjRowIndex();
  int ttm = smps.getTotRows();
  int ttn = smps.getTotCols();

  // count the number of nonzero elements
  int ttnz = smps.countNonzeros(smps.getSmpsTree());

  assert(ttm >= smps.getRows() - 1);
  assert(ttn >= smps.getCols());

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
    int offsetChild;

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

	// stop early before starting a new period
	if (data.rwnmbs[cIndex] >= smps.getBegPeriodRow(period + 1))
      	  break;

	// objective row
	if (data.rwnmbs[cIndex] == iob) {

	  obj[col] = data.acoeff[cIndex] * node->probNode();
	}
	else {
	  acoeff[dIndex] = data.acoeff[cIndex];
	  rwnmbs[dIndex] = data.rwnmbs[cIndex] + offset;

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

      // store the current row index of core
      sIndex = cIndex;

      // store the number of remaining nonzeros in the column
      snzCol = nnzCol;

      // copy the linking blocks in the current period
      for (int block = 0; block < node->nChildren(); ++block) {

	const Node *fChild = node->getChild(0);

	// determine the offset in row numbers for the children nodes
	offsetChild = fChild->firstRow() + block * fChild->nRows()
	  - smps.getBegPeriodRow(period + 1);

	// restore the row index of core and the number of nonzeros
	cIndex = sIndex;
	nnzCol = snzCol;

	// for all the remaining nonzeros in the column
	while (nnzCol > 0) {

	  acoeff[dIndex] = data.acoeff[cIndex];
	  rwnmbs[dIndex] = data.rwnmbs[cIndex] + offsetChild;
#ifdef DEBUG_MATRIX
	  printf("%2d (%2d)| %10f  %d\n",
		 dIndex, cIndex, acoeff[dIndex], rwnmbs[dIndex]);
#endif

	  assert(rwnmbs[dIndex] >= 0);
	  assert(rwnmbs[dIndex] < ttm);

	  ++cIndex;
	  ++dIndex;
	  --nnzCol;
	}
      }

      // go to the next column
      ++col;

    } // end while

  } while (node = node->next());

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

  return PROB;
}

/** Setup the right-hand side */
int setupRhs(ProbData *PROB, const Smps &smps) {

  const Node *node = smps.getRootNode();

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

  }  while (node = node->next());

  PROB->rhs = rhs;
  PROB->rws = rws;

  return 0;
}

/** Apply the scenario corrections */
int applyScenarios(ProbData *PROB, const Smps &smps) {

  int period, firstEntry, lastEntry;
  int row, col, pdr, pdc, index;
  const int *entryRow = smps.getEntryRow();
  const int *entryCol = smps.getEntryCol();
  const double *entryVal = smps.getEntryVal();
  const Node *node = smps.getRootNode(), *scNode;

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
	     scNode, scen + 1, firstEntry, lastEntry - 1);
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

	// if the change affects the matrix
	else if (pdr >= 0 && pdc >= 0 && pdr == period) {

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
	}
#endif
      }

      // walk up the tree to the parent node at the previous period
      scNode = scNode->parent();
      --period;
    }
  } while (node = node->next());

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
  delete PROB;

  return 0;
}

/** Print the solution */
void printSolution(const Node *root,
		   double *primal, double *dual,
		   double *slacks, double *rcosts) {

  int end;
  const Node *node = root;

  // return immediately if there is no root node
  if (!node)
    return;

  // print the solution in breadth-first order
  queue<const Node*> qNodes;

  qNodes.push(node);

  while (!qNodes.empty()) {

    node = qNodes.front();
    qNodes.pop();

    printf("\t---   Node %2d (%dx%d)  ---\n", node->name(),
	   node->nRows(), node->nCols());

    end = node->firstRow() + node->nRows();
    for (int i = node->firstRow(); i < end; ++i) {
      printf("Row %d:  Slack = %10f  Dual = %10f\n", i, slacks[i], dual[i]);
    }

    end = node->firstCol() + node->nCols();
    for (int j = node->firstCol(); j < end; ++j) {
      printf("Column %d:  Value = %10f  Reduced cost = %10f\n",
	     j, primal[j], rcosts[j]);
    }

    // add the children to the queue of nodes to print
    for (int i = 0; i < node->nChildren(); ++i)
      qNodes.push(node->getChild(i));
  }
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
