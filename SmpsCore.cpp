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
  begPeriodCol(NULL),
  nzPeriod(NULL) {
}

/** Destructor */
SmpsCore::~SmpsCore() {

  delete[] rwname;
  delete[] clname;
  delete[] rwstat;
  delete[] varType;
  delete[] hdrwcd;
  delete[] hdclcd;
  delete[] lnkrwcd;
  delete[] lnkclcd;
  delete[] rwnmbs;
  delete[] clpnts;
  delete[] rwhead;
  delete[] links;
  delete[] clnmbs;
  delete[] acoeff;
  delete[] ranges;
  delete[] rhs;
  delete[] bup;
  delete[] blo;
  delete[] begPeriodRow;
  delete[] begPeriodCol;
  delete[] nzPeriod;
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
	++nRows;
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

  rwname  = new char[8 * (maxm + 2)];
  clname  = new char[8 * (maxn + 2)];
  rwstat  = new int[maxm];
  rwhead  = new int[maxm];
  varType = new int[maxn];
  rwnmbs  = new int[maxnza];
  links   = new int[maxnza];
  clnmbs  = new int[maxnza];
  clpnts  = new int[maxn + 1];
  hdrwcd  = new int[maxm + 1];
  hdclcd  = new int[maxn + 1];
  lnkrwcd = new int[maxm + 1];
  lnkclcd = new int[maxn + 1];
  acoeff  = new double[maxnza];
  ranges  = new double[maxm];
  rhs     = new double[maxm];
  bup     = new double[maxn];
  blo     = new double[maxn];
  double *relt = new double[maxn];
  int    *irow = new int[maxn];
  memset(rwname, 0, 8 * (maxm + 2) * sizeof(char));
  memset(clname, 0, 8 * (maxn + 2) * sizeof(char));

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

  // change the numbering convention from fortran to C
  processCore();

 TERMINATE:

  // clean up
  delete[] relt;
  delete[] irow;

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
  for (int i = 0; i < nPeriods; ++i) {

#ifdef DEBUG_TIME_FILE
    cout << "Searching match for row: >" << begPeriodRowName[i] << "<" << endl;
#endif

    // find row-period begin
    for (int j = 0; j < nRows; ++j) {

      if (rowNames[j] == begPeriodRowName[i]) {

	// allow the objective row to be a period start in the time file
	if (j == objRow)
	  ++j;

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
    for (int j = 0; j < nCols; ++j) {

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

  // take care of slacks
  modifyCore();

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
  char ttt[9];

  stream.str(rowname);
  for (int i = 0; i < nRows + 1; ++i) {   // +1 for objective row

    // read exactly 8 characters, at the end of which put the string delimiter
    stream.read(ttt, 8);
    ttt[8] = '\0';

    // create a stream out of the c string to get rid of extra whitespace
    t.str(ttt);
    t >> tmp;
    t.clear();

    // store the cleaned up string in the vector
    rowNames.push_back(tmp);
  }

  stream.clear();
  stream.str(colname);
  for (int i = 0; i < nCols; ++i) {

    // read exactly 8 characters, at the end of which put the string delimiter
    stream.read(ttt, 8);
    ttt[8] = '\0';

    // create a stream out of the c string to get rid of extra whitespace
    t.str(ttt);
    t >> tmp;
    t.clear();

    colNames.push_back(tmp);
  }

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

  for (int i = 0; i < (int) periodNames.size(); ++i) {

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

/**
 *  Determine the period of the given row.
 *
 *  Compares the given row with the period starts stored in begPeriodRow,
 *  and returns the correct period  in [0, nPeriods-1], or -1 if row refers
 *  to the objective row.
 */
int SmpsCore::getRowPeriod(const int row) const {

  int period = -1;
  while (begPeriodRow[period + 1] <= row)
    ++period;

  return period;
}

/**
 *  Determine the period of the given column.
 *
 *  Compares the given column with the period starts stored in begPeriodCol,
 *  and returns the correct period in [0, nPeriods-1], or -1 if col refers
 *  to the right-hand side column.
 */
int SmpsCore::getColPeriod(const int col) const {

  int period = -1;
  while (begPeriodCol[period + 1] <= col)
    ++period;

  return period;
}

/** Set the linked list for row access */
void SmpsCore::setRowsLinkedList() {

  int el, rw;

  // reset the row headers
  for (int row = 0; row < nRows; ++row)
    rwhead[row] = -1;

  // for all columns
  for (int col = 0; col < nCols; ++col) {

    // for all nonzero elements in this column
    for (el = clpnts[col]; el < clpnts[col + 1]; ++el) {
      rw = rwnmbs[el];
      links[el]  = rwhead[rw];
      rwhead[rw] = el;
      clnmbs[el] = col;
    }
  }
}

/** Allocate and return an array for the objective row */
double* SmpsCore::getObjRow() const {

  double *coreObj = new double[nCols];
  memset(coreObj, 0, nCols * sizeof(double));

  int row = rwhead[objRow];
  while (row >= 0) {
    coreObj[clnmbs[row]] = acoeff[row];
    row = links[row];
  }

  return coreObj;
}

/** Retrieve the sparse representation of the matrix */
const SparseData SmpsCore::getSparseData() const {

  SparseData data = {acoeff, rwnmbs, clpnts, clnmbs};

  return data;
}

/**
 *  Change the numbering convention to C.
 *
 *  Changes the core data from Fortran to C numbering convention, orders
 *  elements within columns, and sets up the row linked lists.
 */
void SmpsCore::processCore() {

  int i, el, d1;
  double d2;

  // change to C numbering convention
  for (i = 0; i <= nCols; ++i)
    --clpnts[i];
  for (i = 0; i < nza; ++i)
    --rwnmbs[i];
  --objRow;

  // order entries within columns to have ascending row numbers
  for (int j = 0; j < nCols; ++j) {

    // for all elements in this column
    for (el = clpnts[j]; el < clpnts[j + 1] - 1; ++el) {
      for (i = el + 1; i < clpnts[j + 1]; ++i) {

	if (rwnmbs[el] == rwnmbs[i])
	  printf("Warning: duplicate entry in core matrix: col %d\n", j);

	// check if the elements are not in the correct order
	if (rwnmbs[el] > rwnmbs[i]) {

	  // do a swap
	  d1 = rwnmbs[el];
	  rwnmbs[el] = rwnmbs[i];
	  rwnmbs[i]  = d1;

	  d2 = acoeff[el];
	  acoeff[el] = acoeff[i];
	  acoeff[i]  = d2;
	}
      }
    }
  }

  // set row linked lists for the core matrix
  setRowsLinkedList();
}

/** Count the number of nonzero elements in each period block */
void SmpsCore::countNzPeriodBlocks() {

  int rowPeriod, colPeriod;

  // don't execute it again if it was already run
  if (nzPeriod)
    return;

  // allocate and initialize the vector
  nzPeriod = new int[nPeriods * nPeriods];
  memset(nzPeriod, 0, nPeriods * nPeriods * sizeof(int));

  // for all columns
  for (int i = 0; i < nCols; ++i) {

    // keep track of current column block
    colPeriod = getColPeriod(i) * nPeriods;

    // for all nonzeros in this column
    for (int k = clpnts[i]; k < clpnts[i + 1]; ++k) {

      rowPeriod = getRowPeriod(rwnmbs[k]);
      if (rowPeriod >= 0)
	nzPeriod[rowPeriod + colPeriod]++;
    }
  }

#ifdef DEBUG
  printf("Nonzeros in period blocks of the core matrix:\n");
  for (int i = 0; i < nPeriods; ++i) {
    for (int j = 0; j < nPeriods; ++j)
      printf("  %6d", nzPeriod[i + j * nPeriods]);
    printf("\n");
  }
#endif
}

/**
 *  Convert the inequality constraints to equalities by introducing slacks.
 *
 *  Modify the CORE matrix to convert all ">=" and "<=" constraints into
 *  "=" constraints by introducing slack variables. Slacks are always added
 *  to the last period that is used for this row.
 */
void SmpsCore::modifyCore() {

  int i, j, nSlacks;
  int new_col, new_nnz;

  // count number of inequality rows in CORE matrix
  for (i = 0, nSlacks = 0; i < nRows; ++i)
    if (rwstat[i] == 2 || rwstat[i] == 3) ++nSlacks;

  // new number of columns and nonzeros
  nCols += nSlacks;
  nza   += nSlacks;

  // allocate the resized arrays
  int *new_clpnts = new int[nCols + 1];
  int *new_rwnmbs = new int[nza];
  int *new_stavar = new int[nCols];
  int *new_hdclcd = new int[nCols];
  int *new_lnclcd = new int[nCols];
  int *new_begPeriodCol = new int[nPeriods + 1];
  char *new_clname   = new char[8 * nCols]; // was a calloc
  double *new_acoeff = new double[nza];
  double *new_blo    = new double[nCols];
  double *new_bup    = new double[nCols];
  memset(new_hdclcd, 0, nCols * sizeof(int));

  nSlacks = 0;
  new_col = 0;
  new_nnz = 0;
  new_begPeriodCol[0] = begPeriodCol[0];

  // for all periods
  for (int per = 0; per < nPeriods; ++per) {

    // copy columns of this period
    for (i = begPeriodCol[per]; i < begPeriodCol[per + 1]; ++i) {

      new_clpnts[new_col] = new_nnz;
      new_blo[new_col] = blo[i];
      new_bup[new_col] = bup[i];
      new_stavar[new_col] = varType[i];
      for (j = 0; j < 8; ++j)
	new_clname[8 * new_col + j] = clname[8 * i + j];
      for (j = clpnts[i]; j < clpnts[i + 1]; ++j) {
	new_acoeff[new_nnz] = acoeff[j];
	new_rwnmbs[new_nnz] = rwnmbs[j];
	++new_nnz;
      }
      ++new_col;
    }

    // add slacks for all inequality rows of this period
    for (j = begPeriodRow[per]; j < begPeriodRow[per + 1]; ++j) {

      if (rwstat[j] == 2 || rwstat[j] == 3) {

	// ">=" constraint => add -1 | "<=" constraint => add +1
	new_blo[new_col] = 0.0;
	new_bup[new_col] = 1e31;
	new_clpnts[new_col] = new_nnz;
	new_stavar[new_col] = 1;
	sprintf(&(new_clname[8 * new_col]), "SK%05d", j);
	new_clname[8 * new_col + 7] = ' '; // delete the \0
	new_acoeff[new_nnz] = (rwstat[j] == 2) ? -1.0 : 1.0;
	new_rwnmbs[new_nnz] = j;

	// now consider the row to be an equality
	rwstat[j] = 1;
	++new_col;
	++new_nnz;
	++nSlacks;
      }
    }
    new_begPeriodCol[per + 1] = begPeriodCol[per + 1] + nSlacks;
  }

  // last column definitions
  new_clpnts[new_col] = new_nnz;

  char name[9];
  int kcode;
  int iolog;

  // reset the column names linked list
  for(i = 0; i < nCols; ++i) {

    strncpy(name, new_clname + 8 * i, 8);
    mycode_(&iolog, name, &kcode, &nCols);
    new_lnclcd[i] = new_hdclcd[kcode - 1];
    new_hdclcd[kcode - 1] = i + 1;
  }

  // free the old pointers
  delete[] acoeff;
  delete[] blo;
  delete[] bup;
  delete[] clpnts;
  delete[] rwnmbs;
  delete[] varType;
  delete[] hdclcd;
  delete[] lnkclcd;
  delete[] begPeriodCol;
  delete[] clname;

  // set the new pointers
  acoeff = new_acoeff;
  blo    = new_blo;
  bup    = new_bup;
  clpnts = new_clpnts;
  rwnmbs = new_rwnmbs;
  clname = new_clname;
  hdclcd = new_hdclcd;
  lnkclcd = new_lnclcd;
  varType = new_stavar;
  begPeriodCol = new_begPeriodCol;

#ifdef WITH_MPI
  if (IS_ROOT_PAR)
  printf("Added %d slacks to core matrix.\n", nSlacks);
#endif

#ifdef DEBUG
  printf("Period starts are now:\n"
	 " per    row    col\n");
  for (i = 0; i < nPeriods; ++i)
    printf("%4d %6d %6d\n", i + 1, begPeriodRow[i], begPeriodCol[i]);
  printf(" end %6d %6d\n", begPeriodRow[i], begPeriodCol[i]);
#endif

  // recreate the row linked list for the core matrix
  delete[] clnmbs;
  delete[] links;
  clnmbs = new int[nza + 1];
  links  = new int[nza + 1];
  setRowsLinkedList();
}
