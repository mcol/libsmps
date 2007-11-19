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
#include <sstream>
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
  varType(NULL),
  rwname(NULL),
  clname(NULL),
  hdrwcd(NULL),
  hdclcd(NULL),
  lnkrwcd(NULL),
  lnkclcd(NULL),
  nPeriods(0),
  begPeriodRow(NULL),
  begPeriodCol(NULL) {
}

/** Destructor */
SmpsCore::~SmpsCore() {

  if (rwname)
    free(rwname);
  if (clname)
    free(clname);
  if (rwstat)
    free(rwstat);
  if (varType)
    free(varType);
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
  delete[] begPeriodRow;
  delete[] begPeriodCol;
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

  // check that the file exists
  char core[100] = "";
  strcpy(core, coreFile.c_str());
  if (access(core, R_OK)) {
    cerr << "Could not open file '"  << coreFile << "'." << endl;
    return ERROR_FILE_NOT_FOUND;
  }

  countRows();
  maxm = nRows + 1, maxn = 5 * nRows;
  maxnza = (maxm*maxn/1000 > maxn*10) ? maxm*maxn/1000 : maxn*10;

  rwname  = (char *) calloc(8*(maxm+2), sizeof(char));
  clname  = (char *) calloc(8*(maxn+2), sizeof(char));
  rwstat  = (int *)  calloc(maxm, sizeof(int));
  varType = (int *)  calloc(maxn, sizeof(int));
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
  RDMPS1(&rv, &iqp, &maxm, &maxn, &maxnza, &maxnzq, &nRows, &nCols, &nza, &nzq,
	 &objRow, &iolog, &big, &dlobnd, &dupbnd, &objConstant,
	 namec, nameb, namran, nambnd, nammps, core,
	 qdiag, qclpts, qrwnbs, qcoeff,
	 rwname, clname, varType, rwstat, hdrwcd, lnkrwcd, hdclcd, lnkclcd,
	 rwnmbs, clpnts, irow, acoeff, rhs, ranges, bup, blo, relt);

  if (rv)
    goto TERMINATE;

  // convert the character arrays into vectors of strings
  convertNames(rwname, clname);

 TERMINATE:

  // clean up
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
  vector<string> begPeriodRowName;
  vector<string> begPeriodColName;
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
      begPeriodRowName.push_back(rowName);
      begPeriodColName.push_back(colName);
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

  bool found = false;
  begPeriodRow = new int[nPeriods + 1];
  begPeriodCol = new int[nPeriods + 1];

  // try to match period data
  for (int i = 0; i < nPeriods; i++) {

#ifdef DEBUG_TIME_FILE
    cout << "Searching match for row: >" << begPeriodRowName[i] << "<" << endl;
#endif

    // find row-period begin
    for (int j = 0; j < nRows; j++) {

      if (rowNames[j] == begPeriodRowName[i]) {

	// allow the objective row to be a period start in the time file
	if (j == objRow)
	  j++;

	// found first row of period i
	begPeriodRow[i] = j;
	found = true;
	break;
      }
    }

    if (!found) {
      printf("Not found Row corresponding to start of period %d\n", i);
      rv = 1;
      goto TERMINATE;
    }

#ifdef DEBUG_TIME_FILE
    cout << "Searching match for col: >" << begPeriodColName[i] << "<" << endl;
#endif

    // find col-period begin
    found = false;
    for (int j = 0; j < nCols; j++) {

      if (colNames[j] == begPeriodColName[i]) {

	// found first col of period i
	begPeriodCol[i] = j;
	found = true;
	break;
      }
    }

    if (!found) {
      printf("Not found Col corresponding to start of period %d\n", i);
      rv = 1;
      goto TERMINATE;
    }
  }

  // set the last elements
  begPeriodRow[nPeriods] = nRows;
  begPeriodCol[nPeriods] = nCols;

 TERMINATE:

  return rv;
}

/**
 *  Convert an array of names from rdmps1 into a vector of strings.
 *
 *  Each name is exactly 8 characters long (rdmps1 deals only with
 *  fixed format), so names may be padded with whitespace.
 *  This function extracts each of the names and stores is as an entry
 *  in a vector. During the extraction, the additional whitespace is
 *  removed.
 */
int SmpsCore::convertNames(const char *rowname, const char *colname) {

  string tmp;
  stringstream stream(stringstream::in);
  stringstream t(stringstream::in);
  char *ttt = new char[9];

  stream.str(rowname);
  for (int i = 0; i < nRows + 1; i++) {   // +1 for objective row

    // read exactly 8 characters, at the end of which put the string delimiter
    stream.read(ttt, 8);
    ttt[8] = '\0';

    // create a stream out of the c string to get rid of extra whitespace
    t.str(ttt);
    t >> tmp;
    t.clear();

    // store the clean up string in the vector
    rowNames.push_back(tmp);
  }

  stream.clear();
  stream.str(colname);
  for (int i = 0; i < nCols; i++) {

    // read exactly 8 characters, at the end of which put the string delimiter
    stream.read(ttt, 8);
    ttt[8] = '\0';

    // create a stream out of the c string to get rid of extra whitespace
    t.str(ttt);
    t >> tmp;
    t.clear();

    colNames.push_back(tmp);
  }

  // clean up
  delete[] ttt;

  return 0;
}

/**
 *  Create a single character array from the vector of strings of period names.
 *
 *  The vector of strings of period names read in the time file needs to
 *  be converted into a single character array to be passed to rdstch().
 */
char* SmpsCore::convertPeriodNames() {

  char *periods = new char[8 * nPeriods + 1];
  char *tmp = new char[8 + 1];
  char *pos = periods;

  for (int i = 0; i < (int) periodNames.size(); i++) {

    // extract the character string from the vector
    sprintf(tmp, "%-8s", periodNames[i].c_str());

    // copy the string in the array
    memcpy(pos, tmp, 8 * sizeof(char));
    pos += 8;
  }

  // add the final string delimiter
  periods[8 * nPeriods] = '\0';

  // clean up
  delete[] tmp;

  return periods;
}
