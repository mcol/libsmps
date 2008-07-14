/**
 *  @file SmpsGlpk.h
 *
 *  Definitions for the SMPS interface to Glpk.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
 *
 */

#ifndef _SMPS_GLPK_H_
#define _SMPS_GLPK_H_

#include <glpk.h>
#include "Smps.h"
#include "interface.h"


// Forward declaration
class OptionsGlpk;


/** The SmpsGlpk class */
class SmpsGlpk {

 public:

  /** Constructor */
  SmpsGlpk(string smpsFile);

  /** Destructor */
  ~SmpsGlpk(void);

  /** Read the smps files */
  int read(void);

  /** Solve the problem instance */
  int solve(const OptionsGlpk &opt);

  /** Retrieve the solution information */
  int getSolution(glp_prob *lp, const OptionsGlpk &opt);

 private:

  /** The smps instance to solve */
  Smps smps;

  /** Generate the deterministic equivalent for the smps instance */
  ProbData* generateProblem(void);

};


/** Command line options for the Glpk interface */
class OptionsGlpk : public Options {

  /** Whether an lp file has to be produced */
  int _writeLp;

  /** Whether a barrier solver should be used */
  int _useBarrier;

 public:

  /** Constructor */
  OptionsGlpk(const int argc = 0, const char *argv[] = NULL);

  /** Parse the command line options */
  int parse(void);

  /** Retrieve the value of the writeLp option */
  int writeLp(void) const {
    return _writeLp;
  }

  /** Retrieve the value of the useBarrier option */
  int useBarrier(void) const {
    return _useBarrier;
  }

};

#endif /* _SMPS_GLPK_H_ */
