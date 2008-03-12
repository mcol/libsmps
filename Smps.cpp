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
#include <libgen.h>
#include "Smps.h"


/** Constructor */
Smps::Smps(string smpsFileName) :
  smpsFile(smpsFileName),
  nzPeriod(NULL) {
}

/** Destructor */
Smps::~Smps() {

  delete[] nzPeriod;
}

/** Read the Smps files */
int Smps::read(void) {

  int rv;

  rv = readSmpsFile(smpsFile);
  if (rv)
    return rv;

  rv = readCoreFile();
  if (rv)
    return rv;

  rv = readTimeFile();
  if (rv)
    return rv;

  rv = readStocFile();
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
int Smps::countNonzeros(const Node *rootNode) {

  int nzTotal = 0;
  int nPeriod = getPeriods();

  const Node *node = rootNode;

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
  for (int per = 0, nzPer = 0; per < nPeriod; ++per) {

    // number of nonzeros in the period
    for (int j = 0; j < nPeriod; ++j) {
      nzPer += nzPeriod[per + nPeriod * j];
    }

    // total nonzeros
    nzTotal += nnPer[per] * nzPer;
  }

  // clean up
  delete[] nnPer;

  return nzTotal;
}
