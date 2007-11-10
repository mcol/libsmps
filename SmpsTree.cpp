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
  maxReals(1) {
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
	nBlocks = 0;
      }
    }

    // we reached the ENDATA section
    else if (nTokens == 1 && strcmp(buffer, "ENDATA") == 0) {

      nNodesStage *= (nBlocks + 1);
      maxNodes += nNodesStage;
      maxScens = nNodesStage;

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

      // this block has a different name from the one read before
      if (strcmp(blockName, curBlock) != 0)  {

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

    }

    // we reached the ENDATA section
    else if (nTokens == 1 && strcmp(buffer, "ENDATA") == 0) {

      nNodesStage *= (nBlocks + 1);
      maxNodes += nNodesStage;
      maxScens = nNodesStage;

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
