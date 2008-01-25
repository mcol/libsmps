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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interface.h"

static ProbData* setupMatrix(Smps &smps);
static int setupRhs(ProbData *PROB, Smps &smps);
static int applyScenarios(ProbData *PROB, Smps &smps);


/** Generate the deterministic equivalent problem */
ProbData* setupProblem(Smps &smps) {

  ProbData *PROB;

  printf(" --------------- setupProblem --------------\n");

  // setup the matrix
  PROB = setupMatrix(smps);

  // setup the right-hand side
  setupRhs(PROB, smps);

  // apply the scenario corrections
  applyScenarios(PROB, smps);

  return PROB;
}

/** Setup the constraint matrix */
ProbData *setupMatrix(Smps &smps) {

  int i, node, period, col, block;
  int curCol, nnzCol, snzCol, sIndex;
  int f_cl_nd, f_cl_pd;
  int *rwnmbs, *clpnts;
  double *acoeff, *obj, *blo, *bup;

  int cIndex = 0;   // index in the core matrix
  int dIndex = 0;   // index in the deterministic equivalent

  int nNodes = smps.getNodes();
  int iob = smps.getObjRowIndex();
  int ttm = smps.getFirstRowNode(nNodes);
  int ttn = smps.getFirstColNode(nNodes);

  ProbData *PROB = (ProbData *) calloc(1, sizeof(ProbData));

  // count the number of nonzero elements
  int ttnz = smps.countNonzeros();

  assert(ttm >= smps.getRows() - 1);
  assert(ttn >= smps.getColumns());

  // dimensions of the deterministic equivalent
  printf("Deterministic equivalent matrix is %dx%d, %d nonzeros.\n",
	 ttm, ttn, ttnz);

  // allocate space for the arrays
  obj    = new double[ttn];
  blo    = new double[ttn];
  bup    = new double[ttn];
  acoeff = new double[ttnz];
  rwnmbs = (int *) calloc(ttnz, sizeof(int));
  clpnts = (int *) calloc(ttn + 1, sizeof(int));
  if (obj    == NULL || blo    == NULL || bup    == NULL ||
      acoeff == NULL || rwnmbs == NULL || clpnts == NULL) {
    printf("Memory allocation failed.\n");
    return NULL;
  }
  memset(blo, 0, ttn * sizeof(double));

  SparseData data = smps.getSparseData();

  // for all nodes
  for (node = 0, col = 0; node < smps.getNodes(); ++node) {

    int ordNode = node; //order[node];

    // find which period this node belongs to
    period = smps.getPeriod(ordNode) - 1;

    // first column of the period for the deterministic equivalent and core
    f_cl_nd = smps.getFirstColNode(ordNode);
    f_cl_pd = smps.getBegPeriodCol(period);

    // initialize objective and bounds
    for (i = 0; i < smps.getNColsPeriod(period); ++i) {

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

	  obj[col] = data.acoeff[cIndex] * smps.getProbNode(ordNode);
	}
	else {
	  acoeff[dIndex] = data.acoeff[cIndex];
	  rwnmbs[dIndex] = data.rwnmbs[cIndex] +
	    smps.getFirstRowNode(ordNode) - smps.getBegPeriodRow(period);

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
      for (block = 0; block < smps.getNChildren(ordNode); ++block) {

	// determine the offset in row numbers
	int fChild = smps.getFirstChild(ordNode) - 1;
	int offset = smps.getFirstRowNode(fChild) - smps.getBegPeriodRow(period + 1);

	// restore the row index of core and the number of nonzeros
	cIndex = sIndex;
	nnzCol = snzCol;

	// for all the the remaining nonzeros in the column
	while (nnzCol > 0) {

	  acoeff[dIndex] = data.acoeff[cIndex];
	  rwnmbs[dIndex] = data.rwnmbs[cIndex] +
	    offset + block * smps.getNRowsPeriod(period + 1);
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

  } // end for

  assert(col == ttn);

  // set the last column pointer
  clpnts[col] = dIndex;

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
int setupRhs(ProbData *PROB, Smps &smps) {

  int i, node, period;

  double *rhs = (double *) calloc(PROB->ttm, sizeof(double));
  int    *rws = (int *)    calloc(PROB->ttm, sizeof(int));

  // for all nodes
  for (node = 0; node < smps.getNodes(); ++node) {

    int ordNode = node; //TREE->order[node];

    // find which period the node belongs to
    period = smps.getPeriod(ordNode) - 1;

    // for all rows in this period
    for (i = 0; i < smps.getNRowsPeriod(period); ++i) {

      int index = smps.getFirstRowNode(ordNode) + i;
      assert (index < PROB->ttm);

      rhs[index] = smps.getRhs(smps.getBegPeriodRow(period) + i);
      rws[index] = smps.getRowType(smps.getBegPeriodRow(period) + i);
    }
  }

  PROB->rhs = rhs;
  PROB->rws = rws;

  return 0;
}

/** Apply the scenario corrections */
int applyScenarios(ProbData *PROB, Smps &smps) {

  int scen, node, corr;
  int scNode, period, firstEntry, lastEntry;
  int row, col, pdr, pdc, index;
  const int *entryRow = smps.getEntryRow();
  const int *entryCol = smps.getEntryCol();
  const double *entryVal = smps.getEntryVal();

  // We start from the leaf nodes and traverse the tree up to the root,
  // applying the changes corresponding to the nodes in the scenario
  // only if the correction refers to the same period of the node.

  // find the leaf nodes
  for (node = 0; node < smps.getNodes(); ++node) {

    // skip the non leaf nodes
    if (smps.getNChildren(node) > 0)
      continue;

    scNode = node;
    period = smps.getPeriods() - 1;

    // traverse all nodes from this leaf up to the root
    while (period >= 0) {

      // scenario the current node belongs to
      scen = smps.getScenario(scNode) - 1;

      // index of the first change
      firstEntry = smps.getFirstEntryScen(scen) - 1;
      lastEntry  = firstEntry + smps.getLengthScen(scen);

#ifdef DEBUG_SCEN
      printf("Node %d: scenario %d (entries from %d to %d)\n",
	     scNode, scen + 1, firstEntry, lastEntry - 1);
#endif

      // for all changes affecting this scenario
      for (corr = firstEntry; corr < lastEntry; ++corr) {

	// row and column of core affected by the change
	pdr = smps.getRowPeriod(entryRow[corr] - 1);
	pdc = smps.getColPeriod(entryCol[corr] - 1);

	// if the change affects the objective
	if (pdr < 0 && pdc == period) {

	  col = smps.getFirstColNode(scNode) + entryCol[corr] - 1
	    - smps.getBegPeriodCol(pdc);
	  PROB->obj[col] = entryVal[corr] * smps.getProbNode(scNode);

#ifdef DEBUG_SCEN
	  printf("   Obj entry: core col %3d, det.eq. col %3d, (%f)\n",
		 entryCol[corr], col, entryVal[corr]);
#endif
	}

	// if the change affects the rhs
	else if (pdc < 0 && pdr == period) {

	  row = smps.getFirstRowNode(scNode) + entryRow[corr] - 1
	    - smps.getBegPeriodRow(pdr);
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
	  row = smps.getFirstRowNode(scNode) + entryRow[corr] - 1
	    - smps.getBegPeriodRow(pdr);
	  col = smps.getFirstColNode(scNode) + entryCol[corr] - 1
	    - smps.getBegPeriodCol(pdc);

	  // adjust the column index if the change is in a column that belongs
	  // to the previous period
	  if (pdr != pdc)
	    col += smps.getFirstColNode(smps.getParent(scNode) - 1)
	      - smps.getFirstColNode(scNode);

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
      }

      // walk up the tree to the parent node at the previous period
      scNode = smps.getParent(scNode) - 1;
      --period;
    }
  }

  return 0;
}

/** Free the space allocated for the ProbData structure */
int freeProbData(ProbData *PROB) {

  delete[] PROB->acoeff;
  free(PROB->rwnmbs);
  free(PROB->clpnts);
  delete[] PROB->obj;
  delete[] PROB->blo;
  delete[] PROB->bup;
  free(PROB->rhs);
  free(PROB->rws);
  free(PROB);

  return 0;
}

/** Print the solution */
void printSolution(const NodeInfo *info,
		   double *primal, double *dual,
		   double *slacks, double *rcosts) {

  int i, j, node;

  for (node = 0; node < info->nNodes; ++node) {

    printf("\t---   Node %2d   ---\n", node + 1);

    for (i = info->nRowsNode[node]; i < info->nRowsNode[node + 1]; ++i) {
      printf("Row %d:  Slack = %10f  Dual = %10f\n", i, slacks[i], dual[i]);
    }

    for (j = info->nColsNode[node]; j < info->nColsNode[node + 1]; ++j) {
      printf("Column %d:  Value = %10f  Reduced cost = %10f\n",
	     j, primal[j], rcosts[j]);
    }
  }

#if 0
  int nRows = info->nRowsNode[info->nNodes];
  int nCols = info->nColsNode[info->nNodes];

  for (i = 0, node = 0; i < nRows; ++i) {

    if (i == info->nRowsNode[node]) {
      printf("\t---   Node %2d   ---\n", node + 1);
      ++node;
    }

    printf("Row %d:  Slack = %10f  Dual = %10f\n", i, slacks[i], dual[i]);
  }

  for (j = 0, node = 0; j < nCols; ++j) {

    if (j == info->nColsNode[node]) {
      printf("\t---   Node %2d   ---\n", node + 1);
      ++node;
    }

    printf("Column %d:  Value = %10f  Reduced cost = %10f\n",
	   j, primal[j], rcosts[j]);
  }
#endif
}

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
