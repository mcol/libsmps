/*
 *  SmpsCore.cpp
 *
 *  Structure for the core data of the problem.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include <string.h>
#include <fstream>
#include "Smps.h"
#include "Tokenizer.h"
#include "Utils.h"


/** Constructor */
SmpsCore::SmpsCore(string coreFileName, string timeFileName) :
  coreFile(coreFileName),
  timeFile(timeFileName),
  nRows(0),
  nCols(0),
  nza(0),
  rwnmbs(NULL),
  clpnts(NULL),
  rwhead(NULL),
  links(NULL),
  clnmbs(NULL),
  acoeff(NULL),
  nzq(0),
  qclpts(NULL),
  qrwnbs(NULL),
  qdiag(NULL),
  qcoeff(NULL),
  objRow(0),
  objConstant(0.0),
  blo(NULL),
  bup(NULL),
  rhs(NULL),
  ranges(NULL),
  rwstat(NULL),
  stavar(NULL),
  rwname(NULL),
  clname(NULL),
  hdrwcd(NULL),
  hdclcd(NULL),
  lnkrwcd(NULL),
  lnkclcd(NULL),
  nPeriods(0) {
}

/** Destructor */
SmpsCore::~SmpsCore() {

  if (rwname)
    free(rwname);
  if (clname)
    free(clname);
  if (rwstat)
    free(rwstat);
  if (stavar)
    free(stavar);
  if (hdrwcd)
    free(hdrwcd);
  if (hdclcd)
    free(hdclcd);
  if (lnkrwcd)
    free(lnkrwcd);
  if (lnkclcd)
    free(lnkclcd);
  if (rwnmbs)
    free(rwnmbs);
  if (clpnts)
    free(clpnts);
  if (rwhead)
    free(rwhead);
  if (links)
    free(links);
  if (clnmbs)
    free(clnmbs);
  if (acoeff)
    free(acoeff);
  if (ranges)
    free(ranges);
  if (rhs)
    free(rhs);
  if (bup)
    free(bup);
  if (blo)
    free(blo);
}

/** Count the number of rows declared in the core file */
int SmpsCore::countRows() {

  ifstream core;
  char buffer[SMPS_LINE_MAX];
  char type[SMPS_FIELD_SIZE], name[SMPS_FIELD_SIZE];
  bool foundRows = false;
  int nValuesRead, rv;

  // open the input file
  core.open(coreFile.c_str(), ifstream::in);
  if (core.fail()) {
    cerr << "Could not open file '" << coreFile << "'." << endl;
    return ERROR_FILE_NOT_FOUND;
  }

  // read the file
  while (!core.eof()) {

    // read a line from the file
    rv = readSmpsLine(core, buffer);
    if (rv)
      continue;

    nValuesRead = sscanf(buffer, "%s %s\n", type, name);

    if (nValuesRead == 1) {

      if (strcmp(type, "ROWS") == 0)
	foundRows = true;

      if (strcmp(type, "COLUMNS") == 0)
	break;
    }
    else if (nValuesRead == 2) {
      if (foundRows)
	nRows++;
    }

    else {

    }
  }

  // close the input file
  core.close();

  return 0;
}

/** Read the core file */
int SmpsCore::readCoreFile(string coreFileName) {

  int rv, iqp = 1, iolog = 77;
  double big = 1.e31;
  double dlobnd = 0.0, dupbnd = big;
  int maxm, maxn;
  int maxnza, maxnzq;

  // rdrhs assumes these are initialised like this
  char nameb[10]  = "        ";
  char namec[10]  = "        ";
  char nammps[10] = "        ";
  char nambnd[10] = "        ";
  char namran[10] = "        ";

  // reset SmpsCore::coreFile if a coreFileName has been given
  if (coreFileName != "")
    coreFile = coreFileName;

  countRows();
  maxm = nRows + 1, maxn = 5 * nRows;
  maxnza = (maxm*maxn/1000 > maxn*10) ? maxm*maxn/1000 : maxn*10;

  rwname  = (char *) calloc(8*(maxm+2), sizeof(char));
  clname  = (char *) calloc(8*(maxn+2), sizeof(char));
  rwstat  = (int *)  calloc(maxm, sizeof(int));
  stavar  = (int *)  calloc(maxn, sizeof(int));
  hdrwcd  = (int *)  calloc(maxm+1, sizeof(int));
  hdclcd  = (int *)  calloc(maxn+1, sizeof(int));
  lnkrwcd = (int *)  calloc(maxm+1, sizeof(int));
  lnkclcd = (int *)  calloc(maxn+1, sizeof(int));
  rwnmbs  = (int *)  calloc(maxnza, sizeof(int));
  clpnts  = (int *)  calloc(maxn+1, sizeof(int));
  rwhead  = (int *)  calloc(maxm, sizeof(int));
  links   = (int *)  calloc(maxnza, sizeof(int));
  clnmbs  = (int *)  calloc(maxnza, sizeof(int));
  acoeff  = (double *) calloc(maxnza, sizeof(double));
  ranges  = (double *) calloc(maxm, sizeof(double));
  rhs     = (double *) calloc(maxm, sizeof(double));
  bup     = (double *) calloc(maxn, sizeof(double));
  blo     = (double *) calloc(maxn, sizeof(double));
  double *relt = (double *) calloc(maxn, sizeof(double));
  int    *irow = (int *) calloc(maxn, sizeof(int));

  // read the core file
  char core[100] = "";
  strcpy(core, coreFile.c_str());
  RDMPS1(&rv, &iqp, &maxm, &maxn, &maxnza, &maxnzq, &nRows, &nCols, &nza, &nzq,
	 &objRow, &iolog, &big, &dlobnd, &dupbnd, &objConstant,
	 namec, nameb, namran, nambnd, nammps, core,
	 qdiag, qclpts, qrwnbs, qcoeff,
	 rwname, clname, stavar, rwstat, hdrwcd, lnkrwcd, hdclcd, lnkclcd,
	 rwnmbs, clpnts, irow, acoeff, rhs, ranges, bup, blo, relt);

  free(relt);
  free(irow);

  return rv;
}

/** Read the time file */
int SmpsCore::readTimeFile(string timeFileName) {

  ifstream time;
  char buffer[SMPS_LINE_MAX];
  char rowName[SMPS_FIELD_SIZE], colName[SMPS_FIELD_SIZE];
  char perName[SMPS_FIELD_SIZE];
  bool foundName = false, foundPeriods = false;
  int nTokens, rv = 0;

  // reset SmpsCore::timeFile if a timeFileName has been given
  if (timeFileName != "")
    timeFile = timeFileName;

  // open the time file
  time.open(timeFile.c_str(), ifstream::in);
  if (time.fail()) {
    cerr << "Could not open file '" << timeFile << "'." << endl;
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
      ++nPeriods;
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
