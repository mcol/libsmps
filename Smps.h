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
  SmpsCore(const char *coreFileName = NULL);
  SmpsCore(const char *coreFileName, const char *timeFileName);

  /** Read the time file */
  int readTimeFile(string timeFileName);

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

 private:

  /** Number of rows in the core file */
  int nRows;

  /** Number of periods in the time file */
  int nPeriods;

  /** List of the period names read from the time file */
  vector<string> periodNames;

  /** List of the first row name of each period */
  vector<string> begPeriodRow;

  /** List of the first column name of each period */
  vector<string> begPeriodCol;

  /** Name of the core file */
  char coreFile[SMPS_FILENAME_MAX];

  /** Name of the time file */
  char timeFile[SMPS_FILENAME_MAX];

};


/** The SmpsTree class */
class SmpsTree {

 public:

  /** Constructor */
  SmpsTree(const char *stocFileName = NULL);

  /** Read the stochastic file */
  int readStocFile(const char *stocFileName);

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

  /** Name of the stochastic file */
  char stocFile[SMPS_FILENAME_MAX];

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
class Smps {

 public:

  /** Constructor */
  Smps(const char *smpsFileName = NULL);

  /** Read the Smps files */
  int read(void);

  /** Read the smps input file */
  int readSmpsFile(string smpsFileName);

 private:

  /** Name of the smps input file */
  string smpsFile;

  /** Name of the core file */
  string coreFile;

  /** Name of the time file */
  string timeFile;

  /** Name of the stochastic file */
  string stocFile;

};

#endif /* _SMPS_H_ */