#ifndef _SMPS_H_
#define _SMPS_H_

/** Maximum length of a line */
#define LINE_MAX 200

/** The type of stochastic file */
enum StocType {
  TYPE_INDEP = 1,
  TYPE_BLOCKS,
  LAST_TYPE
};

/** Common error codes */
enum ErrorCodes {
  ERROR_FILE_NOT_FOUND = 1,
  ERROR_WRONG_STOC_TYPE,
  LAST_ERROR
};

/** The maximum between two numbers */
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


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

  /** Maximum number of realisations */
  int maxReals;

  /** Name of the stochastic file */
  char stocFile[100];

  int scanIndepLine(FILE *file);

  int scanBlocksLine(FILE *file);

};

#endif /* _SMPS_H_ */
