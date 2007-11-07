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
Smps::Smps(const char *smpsFileName) :
  nStages(0) {

  if (smpsFileName)
    smpsFile = smpsFileName;
}

/** Read the Smps files */
int Smps::read(void) {

  int rv;

  rv = readSmpsFile();
  if (rv)
    return rv;

  return rv;
}

/** Read the smps input file */
int Smps::readSmpsFile() {

  int rv = 0;
  ifstream smps;

  smps.open(smpsFile.c_str(), ifstream::in);
  if (smps.fail()) {
    cerr << "Could not open file " << smpsFile << endl;
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
