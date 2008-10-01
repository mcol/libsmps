/*
 *  Smps.cpp
 *
 *  Outfacing structure for a problem in Smps format.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <libgen.h>
#include "Smps.h"


/** Constructor */
Smps::Smps(string smpsFileName) :
  smpsFile(smpsFileName),
  Tree(),
  nzPeriod(NULL),
  buildNames(false) {
}

/** Destructor */
Smps::~Smps() {

  delete[] nzPeriod;
}

/**
 *  Read the Smps files.
 *
 *  @param addSlacks
 *         Whether the inequality constraints should be converted into
 *         equalities by adding slacks.
 */
int Smps::read(const bool addSlacks) {

  int rv;

  printf(" --------------- readSmpsFiles -------------\n");

  rv = readSmpsFile(smpsFile);
  if (rv)
    return rv;

  rv = readCoreFile();
  if (rv)
    return rv;

  rv = readTimeFile();
  if (rv)
    return rv;

  // take care of slacks
  if (addSlacks)
    modifyCore();

  rv = readStocFile(Tree);
  if (rv)
    return rv;

  return rv;
}

/** Read the smps input file */
int Smps::readSmpsFile(string smpsFileName) {

  int rv = 0;
  ifstream smps;

  // reset Smps::smpsFile if a smpsFileName has been given
  if (smpsFileName != "")
    smpsFile = smpsFileName;

  // open the input file
  smps.open(smpsFile.c_str(), ifstream::in);
  if (smps.fail()) {
    cerr << "Could not open file '" << smpsFile << "'." << endl;
    return ERROR_FILE_NOT_FOUND;
  }

  // find the path to the problem files
  string path = dirname(const_cast<char *>(smpsFile.c_str()));

  // read the filenames
  smps >> coreFile >> timeFile >> stocFile;

  // insert the path before the filenames
  coreFile = path + "/" + coreFile;
  timeFile = path + "/" + timeFile;
  stocFile = path + "/" + stocFile;

  // close the smps file
  smps.close();

  // check if the files exist and can be read
  if (access(coreFile.c_str(), R_OK)) {
    cerr << "File "  << coreFile << " cannot be read." << endl;
    rv = ERROR_FILE_NOT_FOUND;
  }
  if (access(timeFile.c_str(), R_OK)) {
    cerr << "File "  << timeFile << " cannot be read." << endl;
    rv = ERROR_FILE_NOT_FOUND;
  }
  if (access(stocFile.c_str(), R_OK)) {
    cerr << "File "  << stocFile << " cannot be read." << endl;
    rv = ERROR_FILE_NOT_FOUND;
  }

  return rv;
}

/**
 *  Count the number of nonzeros in the deterministic equivalent matrix.
 *
 *  @note The nonzeros in the objective row are not considered.
 */
int Smps::countNonzeros(const SmpsTree &tree) {

  int nzTotal = 0;
  int nPeriod = getPeriods();

  // get the root node
  const Node *node = tree.getRootNode();

  // leave immediately if there is no root node
  if (!node)
    return 0;

  // count the number of nonzeros in each period block, if not already there
  if (!nzPeriod)
    nzPeriod = countNzPeriodBlocks();

  // number of nodes in each period
  int *nnPer = new int[nPeriod];
  memset(nnPer, 0, nPeriod * sizeof(int));

  // count the number of nodes in each period
  do {

    ++nnPer[node->level()];

  } while (node = node->next());

  // count the number of nonzero elements
  for (int per = 0; per < nPeriod; ++per) {

    // number of nonzeros in the period
    for (int j = 0; j < nPeriod; ++j) {

      nzTotal += nzPeriod[per + nPeriod * j] * nnPer[per];

      // adjust nonzeros for elements above the diagonal
      if (j > per)
	// -1 because it's been added already once at the line above
	nzTotal += nzPeriod[per + nPeriod * j] * (nnPer[j] - 1);
    }
  }

  // clean up
  delete[] nnPer;

  return nzTotal;
}

/**
 *  Set the start rows and columns for each node.
 *
 *  This sets the dimension of each node and the indices of the first
 *  row and column in the deterministic equivalent matrix.
 *
 *  @param tree:
 *         The SmpsTree whose indices need to be updated.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int Smps::setNodeStarts(SmpsTree &tree) {

  int per, rows, cols;
  int ttm = 0, ttn = 0;

  Node *node = tree.getRootNode();

  // leave immediately if there is no root node
  if (!node)
    return 1;

  do {

    per = node->level();
    assert(per < getPeriods());

    // find the dimension of the current node
    rows = cols = 0;
    for (int i = 0; i < node->nLevels(); ++i) {
      rows += getNRowsPeriod(per + i);
      cols += getNColsPeriod(per + i);
    }

    // set the node information
    node->setMatrixPointers(ttm, ttn, rows, cols);

    ttm += rows;
    ttn += cols;

  } while(node = node->next());

  // total number of rows and columns in the deterministic equivalent
  tree.setDimensions(ttm, ttn);

  return 0;
}

/**
 *  Build the row names for the deterministic equivalent.
 *
 *  @return An array of pointers to the row names, which can be NULL
 *          if buildNames is unset.
 *
 *  @note
 *  The name of the objective row is left untouched.
 */
char** Smps::getRowNames(void) const {

  if (!buildNames)
    return NULL;

  const Node *node = getRootNode();
  char **rownames  = new char*[getTotRows()];
  char scname[9];

  // for all nodes in the tree in order
  do {

    sprintf(scname, "_S%03d", node->scenario());
    const int begRowPeriod = getBegPeriodRow(node->level());
    const int firstRowNode = node->firstRow();

    // copy the information for this node
    for (int i = 0; i < node->nRows(); ++i) {

      // build a name for this row
      string rwName = getRowName(begRowPeriod + i) + scname;

      rownames[firstRowNode + i] = new char[20];
      strcpy(rownames[firstRowNode + i], rwName.c_str());
    }

  } while (node = node->next());

  return rownames;
}

/**
 *  Build the column names for the deterministic equivalent.
 *
 *  @return An array of pointers to the column names, which can be NULL
 *          if buildNames is unset.
 */
char** Smps::getColNames(void) const {

  if (!buildNames)
    return NULL;

  const Node *node = getRootNode();
  char **colnames  = new char*[getTotCols()];
  char scname[9];

  // for all nodes in the tree in order
  do {

    sprintf(scname, "_S%03d", node->scenario());
    const int begColPeriod = getBegPeriodCol(node->level());
    const int firstColNode = node->firstCol();

    // copy the information for this node
    for (int i = 0; i < node->nCols(); ++i) {

      // build a name for this column
      string clName = getColName(begColPeriod + i) + scname;

      colnames[firstColNode + i] = new char[20];
      strcpy(colnames[firstColNode + i], clName.c_str());
    }

  } while (node = node->next());

  return colnames;
}
