/**
 *  @file Smps.h
 *
 *  Common declarations for the SMPS classes.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
 *
 */

#ifndef _SMPS_H_
#define _SMPS_H_

#include <fstream>
#include <vector>
#include <string>
#include "Node.h"
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
  TYPE_SCENARIOS,
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


/** Essential sparse matrix representation */
struct SparseData {

  /** Nonzero elements of A */
  const double *acoeff;

  /** Row numbers of A */
  const int *rwnmbs;

  /** Pointers to starts of columns of A */
  const int *clpnts;

  /** Column numbers of elements of A */
  const int *clnmbs;

};


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

  /** Retrieve the number of rows in the core file */
  int getRows(void) const {
    return nRows;
  }

  /** Retrieve the number of columns in the core file */
  int getCols(void) const {
    return nCols;
  }

  /** Retrieve the number of periods in the time file */
  int getPeriods(void) const {
    return nPeriods;
  }

  /** Determine the period for a given row */
  int getRowPeriod(const int row) const;

  /** Determine the period for a given column */
  int getColPeriod(const int col) const;

  /** Return the number of rows of the given period */
  int getNRowsPeriod(const int per) const {
    return begPeriodRow[per + 1] - begPeriodRow[per];
  }

  /** Return the number of columns of the given period */
  int getNColsPeriod(const int per) const {
    return begPeriodCol[per + 1] - begPeriodCol[per];
  }

  /** Return the number of nonzeros in the given period block */
  int getNzPeriod(const int rowBlock, const int colBlock) const {
    return nzPeriod[rowBlock + nPeriods * colBlock];
  }

  /** Retrieve the index of the objective row */
  int getObjRowIndex(void) const {
    return objRow;
  }

  /** Allocate and return an array for the objective row */
  double* getObjRow(void) const;

  /** Return the constraint type */
  int getRowType(const int row) const {
    return rwstat[row];
  }

  /** Return the variable type */
  int getVarType(const int col) const {
    return varType[col];
  }

  /** Return the right-hand side value for the given row */
  double getRhs(const int row) const {
    return rhs[row];
  }

  /** Return the lower bound value for the given variable */
  double getLowerBound(const int col) const {
    return blo[col];
  }

  /** Return the upper bound value for the given variable */
  double getUpperBound(const int col) const {
    return bup[col];
  }

  /** Count the number of nonzeros in each period block */
  void countNzPeriodBlocks(void);

  /** Retrieve the pointer to the beginning of a given row name */
  const char* getBegRowName(const int row) const {
    return &rwname[8 * row];
  }

  /** Retrieve the pointer to the beginning of a given column name */
  const char* getBegColName(const int col) const {
    return &clname[8 * col];
  }

  /** Retrieve the index of the first row of the given period */
  int getBegPeriodRow(const int per) const {
    return begPeriodRow[per];
  }

  /** Retrieve the index of the first column of the given period */
  int getBegPeriodCol(const int per) const {
    return begPeriodCol[per];
  }

  /** Convert the inequality constraints to equalities */
  void modifyCore(void);

  /** Retrieve the sparse representation of the matrix */
  const SparseData getSparseData(void) const;

  /** Give access to SmpsTree to the private members */
  friend class SmpsTree;

 protected:

  /** Name of the core file */
  string coreFile;

  /** Name of the time file */
  string timeFile;

  /** Create a character array from the vector of strings of period names */
  char* convertPeriodNames(void);

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
  int *varType;

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

  /**  Retrieve the index corresponding to the given period name */
  int matchPeriodName(const string& name);

  /** First row index of each period */
  int *begPeriodRow;

  /** First column index of each period */
  int *begPeriodCol;

  /** The number of nonzeros in each period block */
  int *nzPeriod;

  /** List of row names read from the core file */
  vector<string> rowNames;

  /** List of column names read from the core file */
  vector<string> colNames;

  /** Count the number of rows declared in the core file */
  int countRows(void);

  /** Convert the arrays of names from rdmps1 into vectors of strings */
  int convertNames(const char *rowname, const char *colname);

  /** Change the numbering convention to C */
  void processCore(void);

  /** Set the linked list for row access */
  void setRowsLinkedList(void);

  /** Match the names from the time file to those of the core file */
  int findTimeCoreMatches(const vector<string> &begPeriodRowName,
			  const vector<string> &begPeriodColName);

};


/** The SmpsTree class */
class SmpsTree : public SmpsCore {

 public:

  /** Constructor */
  SmpsTree(string stocFileName = "");

  /** Destructor */
  ~SmpsTree();

  /** Read the stochastic file */
  int readStocFile(string stocFileName = "");

  /** Retrieve the root node of the event tree */
  Node* getRootNode(void) const {
    return root;
  }

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

  /** Retrieve the number of rows in the deterministic equivalent */
  int getTotRows(void) const {
    return ttRows;
  }

  /** Retrieve the number of columns in the deterministic equivalent */
  int getTotCols(void) const {
    return ttCols;
  }

  /** Count the number of nonzeros in the deterministic equivalent matrix */
  int countNonzeros(void);

  /** Retrieve the scenario number for the given node */
  int getScenario(const int node) const {
    return scenario[node];
  }

  /** Return the index of the first scenario change for the given scenario */
  const int getFirstEntryScen(const int scen) const {
    return sc_first[scen];
  }

  /** Return the number of changes in the given scenario */
  const int getLengthScen(const int scen) const {
    return sc_len[scen];
  }

  /** Retrieve the pointer to the row indices of the scenario changes */
  const int* getEntryRow(void) const {
    return entryRow;
  }

  /** Retrieve the pointer to the column indices of the scenario changes */
  const int* getEntryCol(void) const {
    return entryCol;
  }

  /** Retrieve the pointer to the values of the scenario changes */
  const double* getEntryVal(void) const {
    return entryVal;
  }

  /** Set the start rows and columns for each node */
  int setNodeStarts(Node *rootNode);

  /** Print the stochastic tree information */
  void printTree(void) const;

 protected:

  /** Name of the stochastic file */
  string stocFile;

 private:

  /** Root node of the event tree */
  Node *root;

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

  /** Number of rows in the deterministic equivalent */
  int ttRows;

  /** Number of columns in the deterministic equivalent */
  int ttCols;

  /** Scenario that node belongs to */
  int *scenario;

  /** Length of scenario correction list */
  int scenLength;

  /** Maximum length of scenario correction list */
  int maxScenLength;

  /** First entry in list for scenario [maxScens] */
  int *sc_first;

  /** Number of entries in list for scenario [maxScens] */
  int *sc_len;

  /* these are in FORTRAN numbering */
  /** Row of CORE affected by change [maxScenLength] */
  int *entryRow;

  /** Column of CORE affected by change [maxScenLength] */
  int *entryCol;

  /** New (changed) entry [maxScenLength] */
  double *entryVal;

  /** Scan a stochastic file in INDEP DISCRETE format */
  int scanIndepType(ifstream &file);

  /** Read a stochastic file in INDEP DISCRETE format */
  int readIndepType(ifstream &file);

  /** Scan a stochastic file in BLOCKS DISCRETE format */
  int scanBlocksType(ifstream &file);

  /** Read a stochastic file in BLOCKS DISCRETE format */
  int readBlocksType(ifstream &file);

  /** Scan a stochastic file in SCENARIOS format */
  int scanScenariosType(ifstream &stoc);

};


/** The Smps class */
class Smps : public SmpsTree {

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
#define RDSTCH      rdstch_
#define MYCODE      mycode_

extern "C" {

/** Read a file in MPS format */
void
RDMPS1(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
       double*, double*, double*, double*,
       char*, char*, char*, char*, char*, const char*,
       double*, int*, int*, double*, char*, char*, int*, int*,
       int*, int*, int*, int*, int*, int*, int*,
       double*, double*, double*, double*, double*, double*);

/** Read a stochastic file */
void
RDSTCH(int*, int*, int*, int*, int*, const char*,
       int*, double*, int*, int*, int*,
       int*, int*, int*, int*, char*, int*, int*, int*, int*,
       int*, int*, double*, int*, int*, char*, char*, int*,
       int*, int*, int*, int*, int*, int*, int*,
       int*, int*, double*, double*, char*,  int*, int*);

/** Compute a hash code */
void
MYCODE(int*, const char*, int*, const int*);

}

#endif /* _SMPS_H_ */
