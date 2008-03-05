/*
 *  SmpsCplex.cpp
 *
 *  SMPS interface to Cplex.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include "SmpsCplex.h"


/** Constructor */
SmpsCplex::SmpsCplex(string smpsFile) :
  smps(smpsFile) {}

/** Destructor */
SmpsCplex::~SmpsCplex() {

}

/** Read the smps files */
int SmpsCplex::read() {

  int rv = smps.read();
  if (rv)
    return rv;

  return rv;
}

/** Generate the deterministic equivalent */
ProbData* SmpsCplex::generateProblem() {

  return setupProblem(smps);
}

/** Constructor */
OptionsCplex::OptionsCplex(const int argc, const char *argv[]) :
  Options(argc, argv),
  _usePresolve(0),
  _useBarrier(0) {
}

/** Parse the command line options */
int OptionsCplex::parse() {

  // add the specialised options
  Options::addOption("-b", "use a barrier solver", &_useBarrier);
  Options::addOption("-p", "turn on the presolve", &_usePresolve);

  // parse the common options
  int rv = Options::parse();

  return rv;
}
