/**
 *  @file SmpsCplex.h
 *
 *  Definitions for the SMPS interface to Cplex.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
 *
 */

#ifndef _SMPS_CPLEX_H_
#define _SMPS_CPLEX_H_

#include <ilcplex/cplex.h>
#include "Smps.h"
#include "interface.h"


// Forward declaration
class OptionsCplex;


/** The SmpsCplex class */
class SmpsCplex {

 public:

  /** Constructor */
  SmpsCplex(string smpsFile);

  /** Destructor */
  ~SmpsCplex(void);

  /** Read the smps files */
  int read(void);

  /** Solve the problem instance */
  int solve(const OptionsCplex &opt);

  /** Retrieve the solution information */
  int getSolution(CPXENVptr env, CPXLPptr lp, const OptionsCplex &opt);

 private:

  /** The smps instance to solve */
  Smps smps;

  /** Generate the deterministic equivalent for the smps instance */
  ProbData* generateProblem(void);

  friend ProbData* setupProblem(Smps &smps);

};


/** Command line options for the Cplex interface */
class OptionsCplex : public Options {

  /** Whether an lp file has to be produced */
  int _writeLp;

  /** Whether the presolver should be used */
  int _usePresolve;

  /** Whether a barrier solver should be used */
  int _useBarrier;

 public:

  /** Constructor */
  OptionsCplex(const int argc = 0, const char *argv[] = NULL);

  /** Parse the command line options */
  int parse(void);

  /** Retrieve the value of the writeLp option */
  int writeLp(void) const {
    return _writeLp;
  }

  /** Retrieve the value of the usePresolve option */
  int usePresolve(void) const {
    return _usePresolve;
  }

  /** Retrieve the value of the useBarrier option */
  int useBarrier(void) const {
    return _useBarrier;
  }

};


#endif /* _SMPS_CPLEX_H_ */
