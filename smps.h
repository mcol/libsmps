#ifndef _SMPS_H_
#define _SMPS_H_

/** Maximum length of a line */
#define LINE_MAX 200

/** Length of each field in the smps files */
#define SMPS_FIELD_SIZE 11

/** The type of stochastic file */
enum StocType {
  TYPE_INDEP = 1,
  TYPE_BLOCKS,
  LAST_TYPE
};

/** Common error codes */
enum ErrorCodes {
  ERROR_FILE_NOT_FOUND = 1,
  ERROR_TIME_FORMAT,
  ERROR_WRONG_STOC_TYPE,
  LAST_ERROR
};

/** The maximum between two numbers */
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


/** The SmpsCore class */
class SmpsCore {

 public:
  SmpsCore(const char *coreFileName);
  SmpsCore(const char *coreFileName, const char *timeFileName);

  /** Count the number of rows declared in the core file */
  int countRows(void);

  /** Retrieve the number of rows in the core file */
  int getRows(void);

  /** Count the number of stages declared in the time file */
  int countStages(void);

  /** Retrieve the number of stages in the problem */
  int getStages(void);

 private:

  /** Number of rows in the core file */
  int nRows;

  /** Number of stages in the problem */
  int nStages;

  /** Name of the core file */
  char coreFile[100];

  /** Name of the time file */
  char timeFile[100];

};


/** The SmpsTree class */
class SmpsTree {

 public:

  /** Constructor */
  SmpsTree(const char *stocFileName);

  /** Parse the stochastic file to retrieve some problem data */
  int getScenarioLength(void);

  /** Retrieve the number of nodes in the event tree */
  int getNodes(void);

  /** Retrieve the maximum number of nodes in the event tree */
  int getMaxNodes(void);

  /** Retrieve the number of scenarios in the event tree */
  int getScens(void);

  /** Retrieve the maximum number of scenarios in the event tree */
  int getMaxScens(void);

  /** Retrieve the number of stages in the event tree */
  int getStages(void);

  /** Retrieve the maximum number of realisations in the event tree */
  int getMaxReals(void);

 private:

  /** Number of nodes in the event tree */
  int nNodes;

  /** Maximum number of nodes in the event tree */
  int maxNodes;

  /** Number of scenarios */
  int nScens;

  /** Maximum number of scenarios */
  int maxScens;

  /** Number of stages in the event tree */
  int nStages;

  /** Maximum number of realisations */
  int maxReals;

  /** Name of the stochastic file */
  char stocFile[100];

  int scanIndepLine(FILE *file);

  int scanBlocksLine(FILE *file);

};

#endif /* _SMPS_H_ */
