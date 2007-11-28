/*
 *  SmpsTree.cpp
 *
 *  Structure for the stochastic data of the problem.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <string.h>
#include <iostream>
#include <fstream>
#include "Smps.h"
#include "Tokenizer.h"
#include "Utils.h"

/**
 *  A BLOCK is a set of entries that will change their values together.
 *  A REALISATION is a set of values for the entries of the BLOCK
 *  together with its probability.
 *
 *  Not all entries are always listed for every BLOCK. The first
 *  time a BLOCK is listed, its list contains all its entries.
 *  After that only entries that are different from the first
 *  listing of this BLOCK are listed again.
 *
 *  So for each BLOCK there are some REALISATIONS and for each
 *  REALISATION there are some ENTRIES.
 *  All realisations of the same block should have the same ENTRIES,
 *  but not always all are listed.
 */

#define DEBUG_SMPSTREE  0

static int
getStocType(char *buffer);

/** Constructor */
SmpsTree::SmpsTree(string stocFileName) :
  stocFile(stocFileName),
  nNodes(0),
  nScens(1),
  nStages(0),
  maxNodes(0),
  maxScens(1),
  maxReals(0),
  parent(NULL),
  nChildren(NULL),
  f_chd(NULL),
  f_rw_nd(NULL),
  f_cl_nd(NULL),
  scenario(NULL),
  period(NULL),
  probnd(NULL),
  scenLength(0),
  maxScenLength(0),
  sc_first(NULL),
  sc_len(NULL),
  entryRow(NULL),
  entryCol(NULL),
  entryVal(NULL) {
}

/** Destructor */
SmpsTree::~SmpsTree() {

  delete[] parent;
  delete[] period;
  delete[] scenario;
  delete[] probnd;
  delete[] nChildren;
  delete[] f_chd;
  delete[] f_rw_nd;
  delete[] f_cl_nd;
  delete[] sc_first;
  delete[] sc_len;
  delete[] entryRow;
  delete[] entryCol;
  delete[] entryVal;
}

/**
 *   Read the stochastic file to retrieve the problem data.
 *
 *   This routine does a quick scan of StocFile and works out the
 *   maximum number of scenarios in the problem, the maximum number
 *   of changes, and the maximum number of nodes in the reduced tree.
 *   After this, a full read of the file is performed, during which the
 *   elements of the SmpsTree structure are set up.
 *
 *   @return 0  Everything is fine.
 *           >0 Error.
 */
int SmpsTree::readStocFile(string stocFileName) {

  ifstream stoc;
  char buffer[SMPS_LINE_MAX];
  char format[SMPS_FIELD_SIZE];
  int foundName = 0;
  int stocType = 0;
  int rv = 0;

  // reset SmpsTree::stocFile if a stocFileName has been given
  if (stocFileName != "")
    stocFile = stocFileName;

  // open the input file
  stoc.open(stocFile.c_str(), ifstream::in);
  if (stoc.fail()) {
    cerr << "Could not open file '" << stocFile << "'." << endl;
    return ERROR_FILE_NOT_FOUND;
  }

  // read the file
  while (!stoc.eof() && !stocType) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

    // find the problem name
    if (!foundName) {
      sscanf(buffer, "%s %*s\n", format);
      if (strcmp(format, "STOCH") == 0)
	foundName = 1;
      continue;
    }

    // determine the type of stochastic file
    stocType = getStocType(buffer);
  }

  // first pass: do a quick scan of the rest of the file
  if (stocType == TYPE_INDEP)
    scanIndepType(stoc);
  else if (stocType == TYPE_BLOCKS)
    scanBlocksType(stoc);
  else
    rv = stocType;

  // close the input file
  stoc.close();

  scenLength = maxScenLength = getMaxReals();

  int *br_sce = new int[maxScens];
  int *iwork1 = new int[maxScenLength];
  int *iwork2 = new int[maxScenLength];
  int *iwork3 = new int[4 * maxScenLength];
  int *iwork4 = new int[maxScenLength];
  char *scenam = new char[8*maxScens];
  double *dwork  = new double[maxScenLength];
  double *prb_rl = new double[maxScenLength];
  char nameb[10];

  parent   = new int[maxNodes];
  period   = new int[maxNodes];
  scenario = new int[maxNodes];
  nChildren= new int[maxNodes];
  f_chd    = new int[maxNodes];
  f_rw_nd  = new int[maxNodes + 1];
  f_cl_nd  = new int[maxNodes + 1];
  probnd   = new double[maxNodes];

  sc_first = new int[maxScens];
  sc_len   = new int[maxScens];
  entryRow = new int[maxScenLength];
  entryCol = new int[maxScenLength];
  entryVal = new double[maxScenLength];

  // initialise to zero
  memset(entryRow, 0, maxScenLength * sizeof(int));
  memset(entryCol, 0, maxScenLength * sizeof(int));
  memset(nChildren, 0, maxNodes * sizeof(int));
  memset(f_chd, 0, maxNodes * sizeof(int));

  int maxRows = nRows;
  int maxCols = nCols;
  int maxPeriods = nPeriods;

  char *perNames = convertPeriodNames();

  // reset SmpsTree::stocFile if a stocFileName has been given
  if (stocFileName != "")
    stocFile = stocFileName;

  // read the stoc file
  char stocfile[100] = "";
  strcpy(stocfile, stocFile.c_str());
  RDSTCH(&rv, &maxScens, &nScens, &maxNodes, &nNodes, stocfile,
	 br_sce, probnd, parent, nChildren, f_chd, scenario,
	 period, &nPeriods, &maxPeriods, perNames,
	 &scenLength, &maxScenLength, entryCol, entryRow,
	 sc_first, sc_len, entryVal,
	 &nCols, &nRows, rwname, clname, &maxCols, &maxRows,
	 hdclcd, hdrwcd, lnkclcd, lnkrwcd, nameb,
	 iwork1, iwork2, iwork3, iwork4, dwork,
	 prb_rl, scenam, begPeriodCol, begPeriodRow);

  // set the start rows and columns for each node
  setNodeStarts();

  /* clean up */
  delete[] perNames;
  delete[] br_sce;
  delete[] iwork1;
  delete[] iwork2;
  delete[] iwork3;
  delete[] iwork4;
  delete[] dwork;
  delete[] prb_rl;
  delete[] scenam;

  return rv;
}

int getStocType(char *buffer) {

  char type[SMPS_FIELD_SIZE], distr[SMPS_FIELD_SIZE];

  sscanf(buffer, "%s %s\n", type, distr);

#if DEBUG_SMPSTREE
  printf(" | Type: %s\n", type);
  printf(" | Distribution: %s\n", distr);
#endif

  // check that the distribution is discrete
  if (strcmp(distr, "DISCRETE") != 0) {
    return TYPE_NOT_IMPLEMENTED;
  }

  if (strcmp(type, "BLOCKS") == 0) {
    return TYPE_BLOCKS;
  }

  else if (strcmp(type, "INDEP") == 0) {
    return TYPE_INDEP;
  }

  else if (strcmp(type, "SCEN") == 0) {
    fprintf(stderr, "Error: SCENARIOS format not implemented.\n");
    return TYPE_NOT_IMPLEMENTED;
  }

  else {
    fprintf(stderr, "Error: Stochastic type not recognised.\n");
    return TYPE_NOT_RECOGNISED;
  }
}

/** Scan a stochastic file in INDEP DISCRETE format */
int SmpsTree::scanIndepType(ifstream &stoc) {

  int rv, nBlocks = 0;
  char buffer[SMPS_LINE_MAX];
  char rowName[SMPS_FIELD_SIZE], curRow[SMPS_FIELD_SIZE] = "";
  char colName[SMPS_FIELD_SIZE], curCol[SMPS_FIELD_SIZE] = "";

  // number of nodes for the current stage
  int nNodesStage = 1;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

    Tokenizer line(buffer);
    int nTokens = line.countTokens();

    if (nTokens == 4 || nTokens == 5) {

      // this is a line of the form
      //  RIGHT     DEMAND1   184.0          (PERIOD2)   0.3
      sscanf(buffer, "%s %s %*f %*s %*f\n", colName, rowName);

      // we have found a new block
      ++nBlocks;

      // check if the name of the block or of the period matches
      if ((strcmp(colName, curCol) != 0) ||
	  (strcmp(rowName, curRow) != 0)) {

	// store the new names of row and period
	strcpy(curRow, rowName);
	strcpy(curCol, colName);

	nNodesStage *= nBlocks;
	maxNodes += nNodesStage;
	++maxReals;
	nBlocks = 0;
      }
    }

    // we reached the ENDATA section
    else if (nTokens == 1 && strncmp(buffer, "ENDATA", 6) == 0) {

      nNodesStage *= (nBlocks + 1);
      maxNodes += nNodesStage;
      maxScens = nNodesStage;
      maxReals = maxScens * maxReals + 1;

      break;
    }

    // we cannot parse this line
    else {
      cerr << "Something wrong with this line?" << endl
	   << ">" << buffer << "<" << endl;
      return ERROR_STOC_FORMAT;
    }
  }

  return 0;
}

/** Read a stochastic file in INDEP DISCRETE format */
int SmpsTree::readIndepType(ifstream &stoc) {

  char buffer[SMPS_LINE_MAX];
  int rv;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;
  }

  return 0;
}

/** Scan a stochastic file in BLOCKS DISCRETE format */
int SmpsTree::scanBlocksType(ifstream &stoc) {

  int rv, nBlocks = 0;
  char buffer[SMPS_LINE_MAX], bl[SMPS_FIELD_SIZE];
  char blockName[SMPS_FIELD_SIZE], curBlock[SMPS_FIELD_SIZE] = "";
  char stageName[SMPS_FIELD_SIZE], curStage[SMPS_FIELD_SIZE] = "";
  bool firstRealBlock = false;

  // number of nodes for the current stage
  int nNodesStage = 1;
  int nRealsBlock = 0;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

    Tokenizer line(buffer);
    int nTokens = line.countTokens();

    if (nTokens == 4) {

      // this is a line of the form
      //  BL BLOCK1    PERIOD2   0.3
      sscanf(buffer, "%s %s %s %*f\n", bl, blockName, stageName);

      // check that the format is the one expected
      if (strcmp(bl, "BL") != 0) {
	cerr << "Something wrong with this line?" << endl
	     << ">" << buffer << "<" << endl;
	return ERROR_STOC_FORMAT;
      }

      // we have found a new block
      ++nBlocks;

      if (firstRealBlock) {
	maxReals += nRealsBlock;
	firstRealBlock = false;
      }

      // this block has a different name from the one read before
      if (strcmp(blockName, curBlock) != 0)  {

	firstRealBlock = true;
	nRealsBlock = 0;

	// store the new name
	strcpy(curBlock, blockName);
	nNodesStage *= nBlocks;
	nBlocks = 0;

	// this block refers to a new period
	if (strcmp(stageName, curStage) != 0) {
	  strcpy(curStage, stageName);
	  maxNodes += nNodesStage;
	}
      }
    }

    // realisation line
    else if (nTokens == 3 || nTokens == 5) {
      ++nRealsBlock;
    }

    // we reached the ENDATA section
    else if (nTokens == 1 && strncmp(buffer, "ENDATA", 6) == 0) {

      nNodesStage *= (nBlocks + 1);
      maxNodes += nNodesStage;
      maxScens = nNodesStage;
      if (firstRealBlock)
	maxReals += nRealsBlock;
      maxReals = maxReals * maxScens + 1;

      break;
    }

    // we cannot parse this line
    else {
      cerr << "Something wrong with this line?" << endl
	   << ">" << buffer << "<" << endl;
      return ERROR_STOC_FORMAT;
    }
  }

  return 0;
}

/** Read a stochastic file in BLOCKS DISCRETE format */
int SmpsTree::readBlocksType(ifstream &stoc) {

  char buffer[SMPS_LINE_MAX];
  int rv;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;
  }

  return 0;
}

/**
 *  Set the start rows and columns for each node.
 *
 *  This sets the arrays f_rw_nd and f_cl_nd with the index in the
 *  deterministic equivalent matrix of the first row and column of
 *  each node. If the nodes have been reordered, then the ordering
 *  array must be provided.
 */
void SmpsTree::setNodeStarts(const int* order) {

  int node, per;

  // number of rows and columns in the deterministic equivalent
  int ttm = 0, ttn = 0;

  // for all nodes
  for (int i = 0; i < nNodes; ++i) {

    if (order)
      node = order[i];
    else
      node = i;
    per = period[node] - 1;

    f_rw_nd[node] = ttm;
    f_cl_nd[node] = ttn;

    ttm += getNRowsPeriod(per);
    ttn += getNColsPeriod(per);
  }

  f_rw_nd[nNodes] = ttm;
  f_cl_nd[nNodes] = ttn;
}

/** Print the stochastic tree information */
void SmpsTree::printTree() const {

  printf("Tree information:\n");
  printf("   node parent scen n_chd  f_chd  per   prob\n");
  for (int i = 0; i < nNodes; i++) {
    printf("  %4i  %4i  %4i  %4i  %4i  %4i   %.4f\n", i + 1,
	   parent[i], scenario[i], nChildren[i],
	   f_chd[i],  period[i],   probnd[i]);
  }
}
