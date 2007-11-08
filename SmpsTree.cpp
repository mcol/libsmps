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
  maxReals(1),
  core(NULL) {
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

  // ensure that we have a core structure attached
  if (!core) {
    printf("No core attached.\n");
    return ERROR_UNATTACHED_CORE;
  }

  // read the SmpsTree::stocFile if no stocFileName has been given
  if (stocFileName == "")
    stocFileName = stocFile;

  // open the input file
  stoc.open(stocFileName.c_str(), ifstream::in);
  if (stoc.fail()) {
    cerr << "Could not open file '" << stocFileName << "'." << endl;
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

  // store the file position at this point
  int filePos = stoc.tellg();

  // first pass: do a quick scan of the rest of the file
  if (stocType == TYPE_INDEP)
    scanIndepType(stoc);
  else if (stocType == TYPE_BLOCKS)
    scanBlocksType(stoc);
  else
    rv = stocType;

  // restore the file position
  stoc.clear();
  stoc.seekg(filePos);

  // second pass: extract all the information from the file
  if (stocType == TYPE_INDEP)
    readIndepType(stoc);
  else if (stocType == TYPE_BLOCKS)
    readBlocksType(stoc);
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

/** Scan a stochastic file in INDEP DISCRETE format */
int SmpsTree::scanIndepType(ifstream &stoc) {

  char buffer[SMPS_LINE_MAX];
  char row[SMPS_FIELD_SIZE], curRow[SMPS_FIELD_SIZE] = "";
  char col[SMPS_FIELD_SIZE], curCol[SMPS_FIELD_SIZE] = "";

  int nChangesBlock = 1;
  int nValuesRead, rv;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

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
	break;
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

  int nBlocks = 0, nRealBlock = 0;
  bool newBlock = false;
  char buffer[SMPS_LINE_MAX];
  char row[SMPS_FIELD_SIZE], curRow[SMPS_FIELD_SIZE] = "";
  char per[SMPS_FIELD_SIZE], curPer[SMPS_FIELD_SIZE] = "";
  char col[SMPS_FIELD_SIZE], rw2[SMPS_FIELD_SIZE];
  int nValuesRead, rv;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

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
	break;
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
