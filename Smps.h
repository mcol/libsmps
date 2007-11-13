/*
 *  Smps.h
 *
 *  Common declarations for the Smps classes.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#ifndef _SMPS_H_
#define _SMPS_H_

#include <fstream>
#include <vector>
#include <string>
using namespace std;


/** Maximum length of a line */
#define SMPS_LINE_MAX 200

/** Maximum length of a filename */
#define SMPS_FILENAME_MAX 100

/** Length of each field in the smps files */
#define SMPS_FIELD_SIZE 11

/** The type of stochastic file */
enum StocType {
  TYPE_INDEP = 1,
  TYPE_BLOCKS,
  TYPE_NOT_IMPLEMENTED,
  TYPE_NOT_RECOGNISED,
  LAST_TYPE
};

/** Common error codes */
enum ErrorCodes {
  ERROR_FILE_NOT_FOUND = 1,
  ERROR_CORE_FORMAT,
  ERROR_TIME_FORMAT,
  ERROR_STOC_FORMAT,
  LAST_ERROR
};

/** The maximum between two numbers */
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


/** The SmpsCore class */
class SmpsCore {

 public:

  /** Constructor */
  SmpsCore(string coreFileName = "", string timeFileName = "");

  /** Destructor */
  ~SmpsCore();

  /** Read the core file */
  int readCoreFile(string coreFileName = "");

  /** Read the time file */
  int readTimeFile(string timeFileName = "");

  /** Count the number of rows declared in the core file */
  int countRows(void);

  /** Retrieve the number of rows in the core file */
  int getRows(void) const {
    return nRows;
  }

  /** Retrieve the number of periods in the time file */
  int getPeriods(void) const {
    return nPeriods;
  }

  /** Give access to SmpsTree to the private members */
  friend class SmpsTree;

 protected:

  /** Name of the core file */
  string coreFile;

  /** Name of the time file */
  string timeFile;

 private:

  /** Number of rows */
  int nRows;

  /** Number of columns */
  int nCols;

  /** Number of nonzeros in A */
  int nza;

  /** Row numbers of A */
  int *rwnmbs;

  /** Pointers to starts of columns of A */
  int *clpnts;

  /** Pointer to last element in row */
  int *rwhead;

  /** Pointer to previous element in row */
  int *links;

  /** Column numbers of elements of A */
  int *clnmbs;

  /** Nonzero elements of A */
  double *acoeff;

  /** Number of nonzeros in Q */
  int nzq;

  /** Column pointers for non-diagonal of Q */
  int *qclpts;

  /** Row indices for non-diag part of Q */
  int *qrwnbs;

  /** Diagonal of Q */
  double *qdiag;

  /** Nonzeros of non-diag part of Q */
  double *qcoeff;

  /** Index of the objective row */
  int objRow;

  /** Constant added to the objective */
  double objConstant;

  /** Variable lower bounds */
  double *blo;

  /** Variable upper bounds */
  double *bup;

  /** Right hand side of problem */
  double *rhs;

  /** Ranges for constraints */
  double *ranges;

  /** Type of constraints */
  int *rwstat;

  /** Type of variables */
  int *stavar;

  /** Array of row names */
  char *rwname;

  /** Array of column names */
  char *clname;

  /** Header to linked list of rows with same code */
  int *hdrwcd;

  /** Header to linked list of columns with same code */
  int *hdclcd;

  /** Linked list of rows with same code */
  int *lnkrwcd;

  /** Linked list of columns with same code */
  int *lnkclcd;

  /** Number of periods in the time file */
  int nPeriods;

  /** List of the period names read from the time file */
  vector<string> periodNames;

  /** First row index of each period */
  int *begPeriodRow;

  /** First column index of each period */
  int *begPeriodCol;

  /** List of row names read from the core file */
  vector<string> rowNames;

  /** List of column names read from the core file */
  vector<string> colNames;

  /** Convert the arrays of names from rdmps1 into vectors of strings */
  int convertNames(const char *rowname, const char *colname);

};


/** The SmpsTree class */
class SmpsTree {

 public:

  /** Constructor */
  SmpsTree(string stocFileName = "");

  /** Read the stochastic file */
  int readStocFile(string stocFileName = "");

  /** Retrieve the number of nodes in the event tree */
  int getNodes(void) const {
    return nNodes;
  }

  /** Retrieve the maximum number of nodes in the event tree */
  int getMaxNodes(void) const {
    return maxNodes;
  }

  /** Retrieve the number of scenarios in the event tree */
  int getScens(void) const {
    return nScens;
  }

  /** Retrieve the maximum number of scenarios in the event tree */
  int getMaxScens(void) const {
    return maxScens;
  }

  /** Retrieve the number of stages in the event tree */
  int getStages(void) const {
    return nStages;
  }

  /** Retrieve the maximum number of realisations in the event tree */
  int getMaxReals(void) const {
    return maxReals;
  }

 protected:

  /** Name of the stochastic file */
  string stocFile;

 private:

  /** Number of nodes in the event tree */
  int nNodes;

  /** Number of scenarios */
  int nScens;

  /** Number of stages in the event tree */
  int nStages;

  /** Maximum number of nodes in the event tree */
  int maxNodes;

  /** Maximum number of scenarios */
  int maxScens;

  /** Maximum number of realisations */
  int maxReals;

  /** Scan a stochastic file in INDEP DISCRETE format */
  int scanIndepType(ifstream &file);

  /** Read a stochastic file in INDEP DISCRETE format */
  int readIndepType(ifstream &file);

  /** Scan a stochastic file in BLOCKS DISCRETE format */
  int scanBlocksType(ifstream &file);

  /** Read a stochastic file in BLOCKS DISCRETE format */
  int readBlocksType(ifstream &file);

};


/** The Smps class */
class Smps : public SmpsCore, public SmpsTree {

 public:

  /** Constructor */
  Smps(string smpsFileName = "");

  /** Read the Smps files */
  int read(void);

  /** Read the smps input file */
  int readSmpsFile(string smpsFileName = "");

 private:

  /** Name of the smps input file */
  string smpsFile;

};

#define RDMPS1      rdmps1_

extern "C" {

void
RDMPS1(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
       double*, double*, double*, double*,
       char*, char*, char*, char*, char*, const char*,
       double*, int*, int*, double*, char*, char*, int*, int*,
       int*, int*, int*, int*, int*, int*, int*,
       double*, double*, double*, double*, double*, double*);
}

#endif /* _SMPS_H_ */
