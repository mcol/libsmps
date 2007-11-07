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
#include "Tokenizer.h"
#include "Utils.h"


/** Constructor */
Smps::Smps(const char smpsFileName[]) :
  nStages(0) {
  // store the file name in the class
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

/** Read the time file */
int Smps::readTimeFile(string timeFileName) {

  ifstream time;
  char buffer[SMPS_LINE_MAX];
  char rowName[SMPS_FIELD_SIZE], colName[SMPS_FIELD_SIZE];
  char perName[SMPS_FIELD_SIZE];
  bool foundName = false, foundPeriods = false;
  int nTokens, rv = 0;

  // open the time file
  time.open(timeFileName.c_str(), ifstream::in);
  if (time.fail()) {
    cerr << "Error: Could not open file " << timeFileName << endl;
    return ERROR_FILE_NOT_FOUND;
  }

  // read the file
  while (!time.eof()) {

    // read a line from the file
    rv = readSmpsLine(time, buffer);
    if (rv)
      continue;

    // find the problem name
    if (!foundName) {
      sscanf(buffer, "%s %*s\n", perName);
      if (strcmp(perName, "TIME") == 0)
	foundName = true;
      continue;
    }

    // find the periods line
    if (!foundPeriods) {
      sscanf(buffer, "%s %*s\n", perName);
      if (strcmp(perName, "PERIODS") == 0)
	foundPeriods = true;
      continue;
    }

    Tokenizer line(buffer);
    nTokens = line.countTokens();

    // normal period declaration line
    if (nTokens == 3) {
      sscanf(buffer, "%s %s %s\n", colName, rowName, perName);
      periodNames.push_back(perName);
      begPeriodRow.push_back(rowName);
      begPeriodCol.push_back(colName);
      ++nStages;
    }

    // we reached the ENDATA section
    else if (nTokens == 1 && strncmp(buffer, "ENDATA", 6) == 0)
      break;

    // we cannot make sense of this line
    else {
      cerr << "Line not recognized (read " << nTokens << " values):\n>"
	   << buffer << "<" << endl;
      continue;
    }
  }

  // we may have reached the end of the file without having found
  // the information we wanted
  if (!foundName || !foundPeriods) {
    cerr << "Problem reading the time file." << endl;
    rv = ERROR_TIME_FORMAT;
  }

  // close the time file
  time.close();

  return rv;
}
