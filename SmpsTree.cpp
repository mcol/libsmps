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
#include <fstream>
#include "smps.h"

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
SmpsTree::SmpsTree(const char *stocFileName) :
  nNodes(0),
  nScens(1),
  nStages(0),
  maxNodes(0),
  maxScens(1),
  maxReals(1) {
  strcpy(stocFile, stocFileName);
}

/** Retrieve the number of nodes in the event tree */
int SmpsTree::getNodes() const {
  return nNodes;
}

/** Retrieve the maximum number of nodes in the event tree */
int SmpsTree::getMaxNodes() const {
  return maxNodes;
}

/** Retrieve the number of scenarios in the event tree */
int SmpsTree::getScens() const {
  return nScens;
}

/** Retrieve the maximum number of scenarios in the event tree */
int SmpsTree::getMaxScens() const {
  return maxScens;
}

/** Retrieve the maximum number of realizations in the event tree */
int SmpsTree::getMaxReals() const {
  return maxReals;
}

/**
 *   Parse the stochastic file to retrieve some problem data.
 *
 *   This routine does a quick scan of StocFile and works out the
 *   maximum number of scenarios in the problem, the maximum number
 *   of changes, and the maximum number of nodes in the reduced tree.
 *
 *   @return 0  Everything is fine.
 *           >0 Error.
 *
 *   @attention
 *   This routine returns an upper bound on the number of nodes.
 */
int SmpsTree::getScenarioLength() {

  ifstream stoc;
  char buffer[SMPS_LINE_MAX];
  char format[SMPS_FIELD_SIZE];
  int foundName = 0;
  int stocType = 0;
  int rv = 0;

  // open the input file
  stoc.open(stocFile, ifstream::in);
  if (stoc.fail()) {
    fprintf(stderr, "Error: Could not open file %s.\n", stocFile);
    return ERROR_FILE_NOT_FOUND;
  }

  // read the file
  while (!stoc.eof() && !stocType) {

    // read a line from the file
    stoc.getline(buffer, SMPS_LINE_MAX);

#if DEBUG_SMPSTREE_BUFFER
    printf(buffer);
#endif

    // skip the asterisk lines
    if (buffer[0] == '*')
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

  // do a quick scan of the rest of the file
  if (stocType == TYPE_INDEP)
    scanIndepLine(stoc);
  else if (stocType == TYPE_BLOCKS)
    scanBlocksLine(stoc);
  else
    rv = stocType;

  // close the input file
  stoc.close();

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

/** Scan the lines of a stochastic file in INDEP DISCRETE format */
int SmpsTree::scanIndepLine(ifstream &stoc) {

  char buffer[SMPS_LINE_MAX];
  char row[SMPS_FIELD_SIZE], curRow[SMPS_FIELD_SIZE];
  char col[SMPS_FIELD_SIZE], curCol[SMPS_FIELD_SIZE];

  int nChangesBlock = 1;
  int nValuesRead;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    stoc.getline(buffer, SMPS_LINE_MAX);

#if DEBUG_SMPSTREE_BUFFER
    printf(buffer);
#endif

    nValuesRead = sscanf(buffer, "%s %s %*f %*s %*f\n", col, row);

    if (nValuesRead == 2) {

      // check if the name of the block or of the period matches
      if ((strcmp(col, curCol) != 0) ||
	  (strcmp(row, curRow) != 0)) {

#if DEBUG_SMPSTREE
	printf(" | Row: %s  Per: %s\n", row, per);
	printf(" | nChangesBlock is: %d\n", nChangesBlock);
#endif

	// update the current values of row and period
	strcpy(curCol, col);
	strcpy(curRow, row);

	maxScens *= MAX(nChangesBlock, 1);
	maxNodes += nChangesBlock * MAX(maxNodes, 1);
	maxReals++;

	nChangesBlock = 1;
      }

      else {
	nChangesBlock += 1;
	//maxReals++;
      }
    }

    else if (nValuesRead == 1) {
      if (strcmp(col, "ENDATA") == 0) {
	maxNodes += nChangesBlock * MAX(maxNodes, 1);
	maxScens *= nChangesBlock;
	maxReals = maxScens * (maxReals - 2) + 5;
      }
    }

    else {
      fprintf(stderr, "Something wrong with this line? (read %d values)\n%s\n",
	      nValuesRead, buffer);
      return ERROR_STOC_FORMAT;
    }
  }

  return 0;
}

/** Scan the lines of a stochastic file in BLOCKS DISCRETE format */
int SmpsTree::scanBlocksLine(ifstream &stoc) {

  int nBlocks = 0, nRealBlock = 0;
  bool newBlock = false;
  char buffer[SMPS_LINE_MAX];
  char row[SMPS_FIELD_SIZE], curRow[SMPS_FIELD_SIZE] = "";
  char per[SMPS_FIELD_SIZE], curPer[SMPS_FIELD_SIZE] = "";
  char col[SMPS_FIELD_SIZE], rw2[SMPS_FIELD_SIZE];
  int nValuesRead;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    stoc.getline(buffer, SMPS_LINE_MAX);

#if DEBUG_SMPSTREE_BUFFER
    printf(buffer);
#endif

    nValuesRead = sscanf(buffer, "%s %s %s %s %*f\n", col, row, per, rw2);

    if (nValuesRead == 4) {

      // beginning of a new block
      if (strcmp(col, "BL") == 0) {
	nBlocks++;

	// special case the first line in each block
	if (newBlock) {
	  maxReals += nRealBlock;
	  newBlock = false;
	}

	// check if the name of the block or of the period matches
	if ((strcmp(row, curRow) != 0) || // row here is the name of the block
	    (strcmp(per, curPer) != 0)) {

#if DEBUG_SMPSTREE
	printf(" | Block: %s  Per: %s\n", curRow, curPer);
	printf(" | Block: %s  Per: %s\n", row, per);
#endif

	  strcpy(curRow, row);
	  strcpy(curPer, per);

	  newBlock = true;
	  nRealBlock = 0;
	  maxScens *= nBlocks;
	  maxNodes += maxScens;
	  nBlocks = 0;
	}
      }

      // normal lines inside the block (with 2 values)
      else {
	nRealBlock++;
      }
    }
    
    // normal lines inside the block (with 1 value)
    else if (nValuesRead == 3) {
      nRealBlock++;
    }
    else if (nValuesRead == 1) {
      if (strcmp(col, "ENDATA") == 0) {
	maxScens = maxScens * (nBlocks + 1);
	maxNodes += maxScens;
	if (newBlock)
	  maxReals += nRealBlock;
	maxReals = maxReals * maxScens;
      }
    }
    else {
      fprintf(stderr, "Something wrong with this line? (read %d values)\n%s\n",
	      nValuesRead, buffer);
      return 1;
    }
  }

  return 0;
}
