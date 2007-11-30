/**
 *  options.c
 *
 *  Command line options.
 *
 *  Andreas Grothey
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include "options.h"
using namespace std;

#ifdef WITH_MPI
#include "oops/parutil.h"
#endif


/** Print a usage message */
static void
showHelpMessage(const char* programName) {

  printf("Usage: %s -f <problem.smps> [options]\n"
	 "   where <problem.smps> is the input SMPS file.\n"
	 "Available options:\n"
	 " -l <number> cutoff level (for multistage programs)\n"
	 " -m   write the deterministic equivalent in MPS format\n"
	 " -s   print the solution\n"
	 " -w   use warmstart\n"
	 " -d   do not solve the problem\n", programName);
}

/** Parse the -d or --dontsolve option */
static int
DontSolveOpt(const int argc, const char *argv[]) {

  for (int arg = 1; arg < argc && argv[arg] != NULL; ++arg) {

    if (strcmp(argv[arg], "-d") == 0 || strcmp(argv[arg], "--dontsolve") == 0)
      return 1;
  }

  return 0;
}

/** Parse the -m or --mps option */
static int
WriteMpsOpt(const int argc, const char *argv[]) {

  for (int arg = 1; arg < argc && argv[arg] != NULL; ++arg) {

    if (strcmp(argv[arg], "-m") == 0 || strcmp(argv[arg], "--mps") == 0)
      return 1;
  }

  return 0;
}

/** Parse the -s or --solution option */
static int
PrintSolutionOpt(const int argc, const char *argv[]) {

  for (int arg = 1; arg < argc && argv[arg] != NULL; ++arg) {

    if (strcmp(argv[arg], "-s") == 0 || strcmp(argv[arg], "--solution") == 0)
      return 1;
  }

  return 0;
}

/** Parse the -o or --output option */
static int
OutputToFileOpt(const int argc, const char *argv[]) {

  for (int arg = 1; arg < argc && argv[arg] != NULL; ++arg) {

    if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0)
      return 1;
  }

  return 0;
}

/** Parse the -w or --warmstart option */
static int UseWarmstartOpt(const int argc, const char *argv[]) {

  for (int arg = 1; arg < argc && argv[arg] != NULL; ++arg) {

    if (strcmp(argv[arg], "-w") == 0 || strcmp(argv[arg], "--warmstart") == 0)
      return 1;
  }

  return 0;
}

/** Get the cutoff level */
static int GetCutoffLevel(const int argc, const char *argv[]) {

  int level = 1;

  for (int arg = 1; arg < argc && argv[arg] != NULL; ++arg) {

    if (strcmp(argv[arg], "-l") == 0 || strcmp(argv[arg], "--level") == 0) {

      if (arg + 1 >= argc) {
	showHelpMessage(argv[0]);
	return -1;
      }
      level = atoi(argv[arg + 1]);
    }
  }

  return level;
}

/**
 *  Get the complete problem description.
 *
 *  The problems must be listed in the following order:
 *   - core file
 *   - time file
 *   - stochastic file
 */
static int
getProblemFiles(const int argc, const char *argv[],
		char **smps, char **core, char **time, char **stoc) {

  for (int arg = 1; arg < argc && argv[arg] != NULL; ++arg) {

    if (strcmp(argv[arg], "-f") == 0) {

      char filename[100], *path;
      FILE *smpsFile;
      int pathlen;

      if (arg + 1 >= argc) {
	showHelpMessage(argv[0]);
	return 1;
      }

      sprintf(filename, "%s", argv[arg + 1]);
      smpsFile = fopen(filename, "r");

      if (!smpsFile) {
	cerr << "Could not open file '" << filename << "'." << endl;
	return 1;
      }

      // allocate space for the filenames
      *smps = (char *) calloc(100, sizeof(char));
      *core = (char *) calloc(100, sizeof(char));
      *time = (char *) calloc(100, sizeof(char));
      *stoc = (char *) calloc(100, sizeof(char));

      // store the smps filename
      sprintf(*smps, "%s", filename);

      // find the path to the problem files
      path = dirname(filename);

      // length of the path string, including the final backslash
      pathlen = strlen(path) + 1;

      // write the initial path
      sprintf(*core, "%s/", path);
      sprintf(*time, "%s/", path);
      sprintf(*stoc, "%s/", path);

      // add the filename
      fscanf(smpsFile, "%s\n", *core + pathlen);
      fscanf(smpsFile, "%s\n", *time + pathlen);
      fscanf(smpsFile, "%s\n", *stoc + pathlen);

      fclose(smpsFile);

      return 0;
    }
  }

  // no -f option given
  fprintf(stderr, "An smps file must be specified with the -f option.\n");

  return 1;
}

/** Parse the command line options */
opt_st* parseOptions(const int argc, const char *argv[]) {

  opt_st *opt = (opt_st *) malloc(sizeof(opt_st));
  int opt_val[6];

  opt->coreFile = NULL;
  opt->timeFile = NULL;
  opt->stocFile = NULL;

#ifdef WITH_MPI
  if (IS_ROOT_PAR) {
#endif

    // the program has been called with not enough options
    if (argc <= 1) {
      showHelpMessage(argv[0]);
      free(opt);
      return NULL;
    }

    int rv = getProblemFiles(argc, argv,
			     &opt->smpsFile, &opt->coreFile,
			     &opt->timeFile, &opt->stocFile);
    if (rv) {
      free(opt);
      return NULL;
    }

    opt_val[0] = GetCutoffLevel(argc, argv);
    opt_val[1] = DontSolveOpt(argc, argv);
    opt_val[2] = PrintSolutionOpt(argc,argv);
    opt_val[3] = OutputToFileOpt(argc, argv);
    opt_val[4] = UseWarmstartOpt(argc, argv);
    opt_val[5] = WriteMpsOpt(argc, argv);

    // check the cutoff level
    if (opt_val[0] < 1) {
      fprintf(stderr, "A positive cutoff level must be specified with -l.\n");
      freeOptions(opt);
      return NULL;
    }

#ifdef WITH_MPI
  }
  MPI_Bcast(opt_val, 6, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  opt->cutoffLevel   = opt_val[0];
  opt->DontSolve     = opt_val[1];
  opt->PrintSolution = opt_val[2];
  opt->OutputToFile  = opt_val[3];
  opt->UseWarmstart  = opt_val[4];
  opt->WriteMps      = opt_val[5];

  return opt;
}

/** Free the space allocated for the options */
void freeOptions(opt_st *opt) {

  if (!opt)
    return;

  free(opt->smpsFile);
  free(opt->coreFile);
  free(opt->timeFile);
  free(opt->stocFile);
  free(opt);
}
