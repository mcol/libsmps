#include <stdio.h>
#include <string.h>
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
#define MPS_FORMAT_SIZE 11

static int
getStocType(char *buffer);

/** Constructor */
SmpsTree::SmpsTree(const char *stocFileName) {
  
  strcpy(stocFile, stocFileName);
  nNodes = 0;
  nScens = 1;
  maxNodes = 0;
  maxScens = 1;
  maxReals = 1;
}

/** Retrieve the number of nodes in the event tree */
int SmpsTree::getNodes() {
  return nNodes;
}

/** Retrieve the maximum number of nodes in the event tree */
int SmpsTree::getMaxNodes() {
  return maxNodes;
}

/** Retrieve the number of scenarios in the event tree */
int SmpsTree::getScens() {
  return nScens;
}

/** Retrieve the maximum number of scenarios in the event tree */
int SmpsTree::getMaxScens() {
  return maxScens;
}

/** Retrieve the maximum number of realizations in the event tree */
int SmpsTree::getMaxReals() {
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
 *
 *   @todo Probably it should not accept 4 parameters, but 2: a class
 *   and the StocFile.
 */
int SmpsTree::getScenarioLength() {

  char buffer[LINE_MAX];
  char format[MPS_FORMAT_SIZE];
  int foundName = 0;
  int stocType = 0;

  // open the input file
  FILE *stoc = fopen(stocFile, "r");
  if (!stoc) {
    fprintf(stderr, "Error: Could not open file %s.\n", stocFile);
    return ERROR_FILE_NOT_FOUND;
  }

  // read all lines of the file
  while (fgets(buffer, LINE_MAX, stoc) != NULL) {

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
    if (!stocType) {
      stocType = getStocType(buffer);

      if (stocType == TYPE_INDEP)
	scanIndepLine(stoc);
      else if (stocType == TYPE_BLOCKS)
	scanBlocksLine(stoc);
      if (stocType < 0)
	return ERROR_WRONG_STOC_TYPE;
    }
  }

  // close the input file
  fclose(stoc);

  return 0;
}

int
getStocType(char *buffer) {

  char type[MPS_FORMAT_SIZE];

  sscanf(buffer, "%s\n", type);

#if DEBUG_SMPSTREE
  printf(" | Type: %s\n", type);
#endif

  if (strcmp(type, "BLOCKS") == 0) {
    return TYPE_BLOCKS;
  }

  else if (strcmp(type, "INDEP") == 0) {
    return TYPE_INDEP;
  }

  else if (strcmp(type, "SCEN") == 0) {
    fprintf(stderr, "Error: SCENARIOS format not implemented.\n");
    return -1;
  }

  else {
    fprintf(stderr, "Error: Stochastic type not recognised.\n");
    return -1;
  }
}

/**
 *  @todo Check that it is DISCRETE.
 */
int SmpsTree::scanIndepLine(FILE *stoc) {

  // check that it is DISCRETE

  char buffer[LINE_MAX];
  char row[MPS_FORMAT_SIZE], curRow[MPS_FORMAT_SIZE];
  char col[MPS_FORMAT_SIZE], curCol[MPS_FORMAT_SIZE];

  int nChangesBlock = 1;
  int nValuesRead;

  // read all lines of the file
  while (fgets(buffer, LINE_MAX, stoc) != NULL) {

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
      return 1;
    }
  }

  return 0;
}

/**
 *  @todo Check that it is DISCRETE.
 */
int SmpsTree::scanBlocksLine(FILE *stoc) {

  int nBlocks = 0, nRealBlock = 0;
  bool newBlock = false;
  char buffer[LINE_MAX];
  char row[MPS_FORMAT_SIZE], curRow[MPS_FORMAT_SIZE] = "";
  char per[MPS_FORMAT_SIZE], curPer[MPS_FORMAT_SIZE] = "";
  char col[MPS_FORMAT_SIZE], rw2[MPS_FORMAT_SIZE];
  int nValuesRead;

  // read all lines of the file
  while (fgets(buffer, LINE_MAX, stoc) != NULL) {

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
