/**
 *  @file options.h
 *
 *  @author Marco Colombo
 *
 */

#ifndef OPTIONS_H
#define OPTIONS_H

/**
 *  Options for the SMPS interface to OOPS.
 */
typedef struct {

  /** Whether the problem should not be solved */
  int DontSolve;

  /** Whether the output has to be redirected to file */
  int OutputToFile;

  /** Cut-off level */
  int cutoffLevel;

  /** Whether a warmstart procedure has to be used */
  int UseWarmstart;

  /** Whether an mps file has to be produced  */
  int WriteMps;

  /** Whether the solution has to be printed to screen */
  int PrintSolution;

  /** Name of the smps file */
  char *smpsFile;

  /** Name of the core file */
  char *coreFile;

  /** Name of the time file */
  char *timeFile;

  /** Name of the stochastic file */
  char *stocFile;

} opt_st;

opt_st* parseOptions(const int argc, const char *argv[]);
void    freeOptions(opt_st *opt);

#endif /* OPTIONS_H */
