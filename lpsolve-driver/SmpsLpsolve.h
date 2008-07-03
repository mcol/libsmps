/**
 *  @file SmpsLpsolve.h
 *
 *  Definitions for the SMPS interface to lpsolve.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
 *
 */

#ifndef _SMPS_LPSOLVE_H_
#define _SMPS_LPSOLVE_H_

#include "Smps.h"
#include "interface.h"
#include "lp_lib.h"


// Forward declaration
class OptionsLpsolve;


/** The SmpsLpsolve class */
class SmpsLpsolve {

 public:

  /** Constructor */
  SmpsLpsolve(string smpsFile);

  /** Destructor */
  ~SmpsLpsolve(void);

  /** Read the smps files */
  int read(void);

  /** Solve the problem instance */
  int solve(const OptionsLpsolve &opt);

  /** Retrieve the solution information */
  int getSolution(lprec *lp, const OptionsLpsolve &opt);

 private:

  /** The smps instance to solve */
  Smps smps;

  /** Generate the deterministic equivalent for the smps instance */
  ProbData* generateProblem(void);

};


/** Command line options for the Lpsolve interface */
class OptionsLpsolve : public Options {

  /** Whether an lp file has to be produced */
  int _writeLp;

  /** Whether the presolver should be used */
  int _usePresolve;

 public:

  /** Constructor */
  OptionsLpsolve(const int argc = 0, const char *argv[] = NULL);

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

};

#endif /* _SMPS_LPSOLVE_H_ */
