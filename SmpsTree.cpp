/*
 *  SmpsTree.cpp
 *
 *  Structure for the stochastic data of the problem.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <string.h>
#include <iostream>
#include <fstream>
#include <queue>
#include "Smps.h"
#include "Tokenizer.h"
#include "Utils.h"

/**
 *  A BLOCK is a set of entries that will change their values together.
 *  A REALISATION is a set of values for the entries of the BLOCK
 *  together with its probability.
 *
 *  Not all entries are always listed for every BLOCK. The first
 *  time a BLOCK is listed, its list contains all its entries.
 *  After that only entries that are different from the first
 *  listing of this BLOCK are listed again.
 *
 *  So for each BLOCK there are some REALISATIONS and for each
 *  REALISATION there are some ENTRIES.
 *  All realisations of the same block should have the same ENTRIES,
 *  but not always all are listed.
 */

static int
getStocType(const char *buffer);

/** Constructor */
SmpsTree::SmpsTree(string stocFileName) :
  stocFile(stocFileName),
  root(NULL),
  nNodes(0),
  nScens(1),
  nStages(0),
  maxNodes(0),
  maxScens(1),
  maxReals(0),
  ttRows(0),
  ttCols(0),
  scenario(NULL),
  scenLength(0),
  maxScenLength(0),
  sc_first(NULL),
  sc_len(NULL),
  entryRow(NULL),
  entryCol(NULL),
  entryVal(NULL) {
}

/** Destructor */
SmpsTree::~SmpsTree() {

  delete root;
  delete[] scenario;
  delete[] sc_first;
  delete[] sc_len;
  delete[] entryRow;
  delete[] entryCol;
  delete[] entryVal;
}

/**
 *   Read the stochastic file to retrieve the problem data.
 *
 *   This routine does a quick scan of StocFile and works out the
 *   maximum number of scenarios in the problem, the maximum number
 *   of changes, and the maximum number of nodes in the reduced tree.
 *   After this, a full read of the file is performed, during which the
 *   elements of the SmpsTree structure are set up.
 *
 *   @return 0  Everything is fine.
 *           >0 Error.
 */
int SmpsTree::readStocFile(string stocFileName) {

  ifstream stoc;
  char buffer[SMPS_LINE_MAX];
  char format[SMPS_FIELD_SIZE];
  int foundName = 0;
  int stocType = 0;
  int rv = 0;

  // index in the array representation of the event tree
  int index = 0;

  // name of the node for the node representation of the event tree
  int nodeName = 1;

  // queue of nodes to be processed
  queue<Node*> qNodes;

  // reset SmpsTree::stocFile if a stocFileName has been given
  if (stocFileName != "")
    stocFile = stocFileName;

  // open the input file
  stoc.open(stocFile.c_str(), ifstream::in);
  if (stoc.fail()) {
    cerr << "Could not open file '" << stocFile << "'." << endl;
    return ERROR_FILE_NOT_FOUND;
  }

  // read the file
  while (!stoc.eof() && !stocType) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

    // find the problem name
    if (!foundName) {
      sscanf(buffer, "%s %*s\n", format);
      if (strcmp(format, "STOCH") == 0)
	foundName = 1;
      continue;
    }

    // determine the type of stochastic file
    stocType = getStocType(buffer);
  }

  // first pass: do a quick scan of the rest of the file
  switch(stocType) {

  case TYPE_INDEP:
    rv = scanIndepType(stoc);
    break;
  case TYPE_BLOCKS:
    rv = scanBlocksType(stoc);
    break;
  case TYPE_SCENARIOS:
    rv = scanScenariosType(stoc);
    break;

  default:
    rv = stocType;
  }

  // close the input file
  stoc.close();

  // early return if something has gone wrong
  if (rv)
    return rv;

  scenLength = maxScenLength = getMaxReals();
  ++maxScenLength;

  int *parent = new int[maxNodes];
  int *f_chd  = new int[maxNodes];
  int *n_chd  = new int[maxNodes];
  int *period = new int[maxNodes];
  int *br_sce = new int[maxScens];
  int *iwork1 = new int[maxScenLength];
  int *iwork2 = new int[maxScenLength];
  int *iwork3 = new int[4 * maxScenLength];
  int *iwork4 = new int[maxScenLength];
  char *scenam = new char[8*maxScens];
  double *probnd = new double[maxNodes];
  double *dwork  = new double[maxScenLength];
  double *prb_rl = new double[maxScenLength];

  scenario = new int[maxNodes];
  sc_first = new int[maxScens];
  sc_len   = new int[maxScens];
  entryRow = new int[maxScenLength];
  entryCol = new int[maxScenLength];
  entryVal = new double[maxScenLength];

  // initialise to zero
  memset(entryRow, 0, maxScenLength * sizeof(int));
  memset(entryCol, 0, maxScenLength * sizeof(int));
  memset(n_chd, 0, maxNodes * sizeof(int));
  memset(f_chd, 0, maxNodes * sizeof(int));

  int maxRows = getRows();
  int maxCols = getCols();
  int maxPers = getPeriods();

  char *perNames = convertPeriodNames();

  // reset SmpsTree::stocFile if a stocFileName has been given
  if (stocFileName != "")
    stocFile = stocFileName;

  // read the stoc file
  char stocfile[100] = "";
  strcpy(stocfile, stocFile.c_str());
  RDSTCH(&rv, &maxScens, &nScens, &maxNodes, &nNodes, stocfile,
	 br_sce, probnd, parent, n_chd, f_chd, scenario,
	 period, &maxPers, &maxPers, perNames,
	 &scenLength, &maxScenLength, entryCol, entryRow,
	 sc_first, sc_len, entryVal,
	 &maxCols, &maxRows, rwname, clname, &maxCols, &maxRows,
	 hdclcd, hdrwcd, lnkclcd, lnkrwcd,
	 iwork1, iwork2, iwork3, iwork4, dwork,
	 prb_rl, scenam, begPeriodCol, begPeriodRow);

  if (rv)
    goto TERMINATE;

  // Convert the event tree from an array-based representation into a
  // node-based one. We go through the nodes in the array and set up a
  // tree in breadth-first order.
  // Note that the scenario information is still kept in arrays.

  // allocate the root node
  root = new Node(nodeName);

  qNodes.push(root);

  while (!qNodes.empty()) {

    // take the first element in the queue
    Node *node = qNodes.front();
    qNodes.pop();

    // add the information to this node
    node->setScen(scenario[index] - 1);
    node->setProb(probnd[index]);

    // allocate its children
    for (int j = 0; j < n_chd[index]; ++j) {

      Node *child = new Node(++nodeName);
      node->addChildNode(child);
      qNodes.push(child);
    }

    // set an initial breadth-first ordering
    node->setNext(qNodes.empty() ? NULL : qNodes.front());

    // consider the next node in the array
    ++index;
  }

  // set the start rows and columns for each node
  rv = setNodeStarts(root);
  if (rv)
    goto TERMINATE;

#ifdef DEBUG
  if (nNodes <= 100)
    printTree(root);
#endif

 TERMINATE:

  // clean up
  delete[] parent;
  delete[] f_chd;
  delete[] n_chd;
  delete[] period;
  delete[] perNames;
  delete[] br_sce;
  delete[] iwork1;
  delete[] iwork2;
  delete[] iwork3;
  delete[] iwork4;
  delete[] dwork;
  delete[] probnd;
  delete[] prb_rl;
  delete[] scenam;

  return rv;
}

/** Find out the format of the stochastic file */
int getStocType(const char *buffer) {

  char type[SMPS_FIELD_SIZE], distr[SMPS_FIELD_SIZE];

  int nValuesRead = sscanf(buffer, "%s %s\n", type, distr);

#ifdef DEBUG_SMPS_TREE
  printf(" | Type: %s\n", type);
#endif

  if (nValuesRead > 1) {
    // check that the distribution is discrete
    if (strcmp(distr, "DISCRETE") != 0 || strcmp(distr, "") == 0 ) {
      fprintf(stderr, "Error: Distribution '%s' not recognised.\n", distr);
      return TYPE_NOT_RECOGNISED;
    }
  }

  if (strcmp(type, "INDEP") == 0) {
    return TYPE_INDEP;
  }

  else if (strcmp(type, "BLOCKS") == 0) {
    return TYPE_BLOCKS;
  }

  else if (strcmp(type, "SCENARIOS") == 0 ||
	   strcmp(type, "SCEN") == 0) {
    return TYPE_SCENARIOS;
  }

  else {
    fprintf(stderr, "Error: Stochastic type not recognised.\n");
    return TYPE_NOT_RECOGNISED;
  }
}

/** Scan a stochastic file in INDEP DISCRETE format */
int SmpsTree::scanIndepType(ifstream &stoc) {

  int rv, nBlocks = 0;
  char buffer[SMPS_LINE_MAX];
  char rowName[SMPS_FIELD_SIZE], colName[SMPS_FIELD_SIZE];
  string curRow = "", curCol = "";

  // number of nodes for the current stage
  int nNodesStage = 1;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

    Tokenizer line(buffer);
    int nTokens = line.countTokens();

    if (nTokens == 4 || nTokens == 5) {

      // this is a line of the form
      //  RIGHT     DEMAND1   184.0          (PERIOD2)   0.3
      sscanf(buffer, "%s %s %*f %*s %*f\n", colName, rowName);

      // we have found a new block
      ++nBlocks;

      // check if the name of the block or of the period matches
      if ((colName != curCol) || (rowName != curRow)) {

	// store the new names of row and period
	curRow = rowName;
	curCol = colName;

	nNodesStage *= nBlocks;
	maxNodes += nNodesStage;
	++maxReals;
	nBlocks = 0;
      }
    }

    // we reached the ENDATA section
    else if (nTokens == 1 && strncmp(buffer, "ENDATA", 6) == 0) {

      nNodesStage *= (nBlocks + 1);
      maxNodes += nNodesStage;
      maxScens = nNodesStage;
      maxReals = maxScens * maxReals + 1;

      break;
    }

    // we cannot parse this line
    else {
      cerr << "Something wrong with this line?" << endl
	   << ">" << buffer << "<" << endl;
      return ERROR_STOC_FORMAT;
    }
  }

  return 0;
}

/** Read a stochastic file in INDEP DISCRETE format */
int SmpsTree::readIndepType(ifstream &stoc) {

  char buffer[SMPS_LINE_MAX];
  int rv;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;
  }

  return 0;
}

/** Scan a stochastic file in BLOCKS DISCRETE format */
int SmpsTree::scanBlocksType(ifstream &stoc) {

  int rv, nBlocks = 0;
  char buffer[SMPS_LINE_MAX], bl[SMPS_FIELD_SIZE];
  char blockName[SMPS_FIELD_SIZE], stageName[SMPS_FIELD_SIZE];
  string curBlock = "", curStage = "";
  bool firstRealBlock = false;

  // number of nodes for the current stage
  int nNodesStage = 1;
  int nRealsBlock = 0;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

    Tokenizer line(buffer);
    int nTokens = line.countTokens();

    if (nTokens == 4) {

      // this is a line of the form
      //  BL BLOCK1    PERIOD2   0.3
      sscanf(buffer, "%s %s %s %*f\n", bl, blockName, stageName);

      // check that the format is the one expected
      if (strcmp(bl, "BL") != 0) {
	cerr << "Something wrong with this line?" << endl
	     << ">" << buffer << "<" << endl;
	return ERROR_STOC_FORMAT;
      }

      // we have found a new block
      ++nBlocks;

      if (firstRealBlock) {
	maxReals += nRealsBlock;
	firstRealBlock = false;
      }

      // this block has a different name from the one read before
      if (blockName != curBlock)  {

	firstRealBlock = true;
	nRealsBlock = 0;

	// store the new name
	curBlock = blockName;
	nNodesStage *= nBlocks;
	nBlocks = 0;

	// this block refers to a new period
	if (stageName != curStage) {
	  curStage = stageName;
	  maxNodes += nNodesStage;
	}
      }
    }

    // realisation line
    else if (nTokens == 3 || nTokens == 5) {
      ++nRealsBlock;
    }

    // we reached the ENDATA section
    else if (nTokens == 1 && strncmp(buffer, "ENDATA", 6) == 0) {

      nNodesStage *= (nBlocks + 1);
      maxNodes += nNodesStage;
      maxScens = nNodesStage;
      if (firstRealBlock)
	maxReals += nRealsBlock;
      maxReals = maxReals * maxScens + 1;

      break;
    }

    // we cannot parse this line
    else {
      cerr << "Something wrong with this line?" << endl
	   << ">" << buffer << "<" << endl;
      return ERROR_STOC_FORMAT;
    }
  }

  return 0;
}

/** Read a stochastic file in BLOCKS DISCRETE format */
int SmpsTree::readBlocksType(ifstream &stoc) {

  char buffer[SMPS_LINE_MAX];
  int rv;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;
  }

  return 0;
}

/** Scan a stochastic file in SCENARIOS format */
int SmpsTree::scanScenariosType(ifstream &stoc) {

  int rv, branchPeriod;
  char buffer[SMPS_LINE_MAX], sc[SMPS_FIELD_SIZE], perName[SMPS_FIELD_SIZE];
  char scenName[SMPS_FIELD_SIZE], fromName[SMPS_FIELD_SIZE];
  maxScens = 0;

  // read the file
  while (!stoc.eof()) {

    // read a line from the file
    rv = readSmpsLine(stoc, buffer);
    if (rv)
      continue;

    Tokenizer line(buffer);
    int nTokens = line.countTokens();

    // scenario declaration line
    if (nTokens == 5) {

      // this is a line of the form
      //  SC SCEN_A    'ROOT'    0.09           PERIOD1
      sscanf(buffer, "%s %s %s %*f %s\n", sc, scenName, fromName, perName);

      // check that the format is the one expected
      if (strcmp(sc, "SC") != 0) {
	cerr << "Something wrong with this line?" << endl
	     << ">" << buffer << "<" << endl;
	return ERROR_STOC_FORMAT;
      }

      // we have found a new scenario
      ++maxScens;

      // find out at what period this scenario branches
      branchPeriod = matchPeriodName(perName);
      if (branchPeriod < 0) {
	cerr << "Period >" << perName << "< not declared in the time file."
	     << endl;
	return ERROR_STOC_FORMAT;
      }

      maxNodes += nPeriods - branchPeriod;
    }

    // realisation line
    else if (nTokens == 3) {

      ++maxReals;
    }

    // we reached the ENDATA section
    else if (nTokens == 1 && strncmp(buffer, "ENDATA", 6) == 0) {

      break;
    }

    // we cannot parse this line
    else {
      cerr << "Something wrong with this line?" << endl
	   << ">" << buffer << "<" << endl;
      return ERROR_STOC_FORMAT;
    }
  }

  return 0;
}

/**
 *  Count the number of nonzeros in the deterministic equivalent matrix.
 *
 *  @note The nonzeros in the objective row are not considered.
 */
int SmpsTree::countNonzeros(const Node *rootNode) {

  int nzTotal = 0;
  int nPeriod = getPeriods();

  const Node *node = rootNode;

  // leave immediately if there is no root node
  if (!node)
    return 1;

  // count the number of nonzeros in each period block, if not already there
  if (!nzPeriod)
    countNzPeriodBlocks();

  // number of nodes in each period
  int *nnPer = new int[nPeriod];
  memset(nnPer, 0, nPeriod * sizeof(int));

  // count the number of nodes in each period
  do {

    ++nnPer[node->level()];

  } while (node = node->next());

  // count the number of nonzero elements
  for (int per = 0; per < nPeriod; ++per) {

    int nzPer = 0;

    // number of nonzeros in the period
    for (int j = 0; j < nPeriod; ++j) {
      nzPer += nzPeriod[per + j * nPeriod];
    }

    // total nonzeros
    nzTotal += nnPer[per] * nzPer;
  }

  // clean up
  delete[] nnPer;

  return nzTotal;
}

/**
 *  Set the start rows and columns for each node.
 *
 *  This sets the dimension of each node and the indices of the first
 *  row and column in the deterministic equivalent matrix.
 *
 *  @param rootNode:
 *         The root node of the tree whose indices need to be updated.
 *  @return 1 If something goes wrong; 0 otherwise.
 */
int SmpsTree::setNodeStarts(Node *rootNode) {

  int per, rows, cols;
  int ttm = 0, ttn = 0;

  Node *node = rootNode;

  // leave immediately if there is no root node
  if (!node)
    return 1;

  do {

    per  = node->level();
    rows = getNRowsPeriod(per);
    cols = getNColsPeriod(per);

    // set the node information
    node->setMatrixPointers(ttm, ttn, rows, cols);

    ttm += rows;
    ttn += cols;

  } while(node = node->next());

  // total number of rows and columns in the deterministic equivalent
  ttRows = ttm;
  ttCols = ttn;

  return 0;
}

/** Print the stochastic tree information */
void SmpsTree::printTree(const Node *rootNode) const {

  const Node *node = rootNode;

  // leave immediately if there is no root node
  if (!node)
    return;

  // queue of nodes to be printed out in breadth-first order
  queue<const Node*> qNodes;

  // start from the root
  qNodes.push(node);

  printf("Tree information:\n");
  printf("   node parent scen  chdn   per   prob   |  rows  cols\n");

  // go through all nodes in the tree in breadth-first order
  while (!qNodes.empty()) {

    // take the first element in the queue
    node = qNodes.front();
    qNodes.pop();

    // print it
    node->print();

    // add its children to the queue
    for (int i = 0; i < node->nChildren(); ++i)
      qNodes.push(node->getChild(i));
  }
}
