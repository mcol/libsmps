/**
 *  @file SmpsHopdm.h
 *
 *  Definitions for the SMPS interface to Hopdm.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
 *
 */

#ifndef _SMPS_HOPDM_H_
#define _SMPS_HOPDM_H_

#define INTS long

#include "Smps.h"
#include "interface.h"
extern "C" {
#include "lp.h"
#include "mps.h"
}

typedef Lp_save Ws;

// Forward declaration
class OptionsHopdm;


/** The SmpsHopdm class */
class SmpsHopdm {

 public:

  /** Constructor */
  SmpsHopdm(string smpsFile);

  /** Read the smps files */
  int read(void);

  /** Solve the problem instance */
  int solve(const OptionsHopdm &opt);

  /** Solve a reduced problem */
  int solveReduced(const OptionsHopdm &opt);

  /** Retrieve the solution information */
  int getSolution(Lp *lp, const OptionsHopdm &opt);

 private:

  /** The smps instance to solve */
  Smps smps;

  /** The reduced event tree */
  SmpsTree rTree;

  /** Generate the deterministic equivalent for the smps instance */
  ProbData* generateProblem(const Node *root);

  /** Store the solution from the reduced problem */
  //  int storeSolution(const PDProblem *pdProb, const SmpsReturn *Ret);

  /** Set up the warmstart point from a reduced-tree solution */
  //  int setupWarmStart(const SmpsReturn *Pb);

};

/*
int hopdmSolve(ProbData *prob, Ws **iterate, const int mode);
void setupWarmStart(const SmpsTree *cTree, const SmpsTree *rTree,
                    Ws **wsIterate, const int *chosen);
*/

/** Command line options for the Hopdm interface */
class OptionsHopdm : public Options {

  /** Whether scenario reduction should be used */
  int _useReduction;

  /** Whether stage aggregation should be performed */
  int _useAggregation;

 public:

  /** Constructor */
  OptionsHopdm(const int argc = 0, const char *argv[] = NULL);

  /** Parse the command line options */
  int parse(void);

  /** Determine whether a warmstart strategy has to be employed */
  int useWarmstart(void) const {
    return _useReduction || _useAggregation;
  }

  /** Retrieve the value of the useReduction option */
  int useReduction(void) const {
    return _useReduction;
  }

  /** Retrieve the value of the useAggregation option */
  int useAggregation(void) const {
    return _useAggregation;
  }

};

#endif /* _SMPS_HOPDM_H_ */
