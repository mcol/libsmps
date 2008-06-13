/**
 *  @file interface.h
 *
 *  Main header for the SMPS interface.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
 *
 */

#ifndef _SMPS_INTERFACE_H
#define _SMPS_INTERFACE_H

#include "Smps.h"
#include "options.h"


/** Data for the deterministic equivalent in sparse format */
struct ProbData {

  /** Number of rows */
  int ttm;

  /** Number of columns */
  int ttn;

  /** Number of nonzeros */
  int ttnz;

  /** Index of the objective row */
  int irobj;

  /** nonzero elements of A */
  double *acoeff;

  /** row numbers of A */
  int *rwnmbs;

  /** pointers to starts of columns of A */
  int *clpnts;

  /** Objective function */
  double *obj;

  /** Variable lower bounds */
  double *blo;

  /** Variable upper bounds */
  double *bup;

  /** Right-hand side vector */
  double *rhs;

  /** Type of constraint */
  int *rws;

  /** Row names */
  char **rwnames;

  /** Column names */
  char **clnames;

};


/** Free the space allocated for a ProbData structure */
int freeProbData(ProbData *data);

/** Setup the deterministic equivalent */
ProbData* setupProblem(Smps &smps);

/** Print the solution */
void printSolution(const Node *root,
		   double *primal, double *dual,
		   double *slacks, double *rcosts);

/** Write the deterministic equivalent in MPS format */
int writeMps(char *filename, ProbData *PROB);

#endif /* _SMPS_INTERFACE_H */
