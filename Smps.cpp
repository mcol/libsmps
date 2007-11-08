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
  smpsFile(smpsFileName) {
}

/** Read the Smps files */
int Smps::read(void) {

  int rv;

  rv = readSmpsFile(smpsFile);
  if (rv)
    return rv;

  return rv;
}

/** Read the smps input file */
int Smps::readSmpsFile(string smpsFileName) {

  int rv = 0;
  ifstream smps;

  smps.open(smpsFileName.c_str(), ifstream::in);
  if (smps.fail()) {
    cerr << "Could not open file " << smpsFileName << endl;
    return ERROR_FILE_NOT_FOUND;
  }

  // find the path to the problem files
  string path = dirname(const_cast<char *>(smpsFileName.c_str()));

  // read the filenames
  smps >> coreFile >> timeFile >> stocFile;

  // insert the path before the filenames
  coreFile = path + "/" + coreFile;
  timeFile = path + "/" + timeFile;
  stocFile = path + "/" + stocFile;

  // set the filenames of the base classes
  SmpsCore::coreFile = coreFile;
  SmpsCore::timeFile = timeFile;
  SmpsTree::stocFile = stocFile;

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
