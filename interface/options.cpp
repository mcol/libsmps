/*
 *  options.cpp
 *
 *  Command line options.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include "options.h"


/** Constructor */
Options::Options(const int argc, const char *argv[]) :
  nArgs(argc),
  vArgs(argv),
  _dontSolve(0),
  _outputToFile(0),
  _writeMps(0),
  _printSolution(0),
  _smpsFile("") {

  // add the common options
  addOption("-s", "print the solution",
	    &_printSolution);
  addOption("-d", "do not solve the problem",
	    &_dontSolve);
  addOption("-m", "write the deterministic equivalent in MPS format",
	    &_writeMps);
}

/** Destructor */
Options::~Options() {
}

/** Parse the command line options */
int Options::parse() {

  const char *programName = vArgs[0];

  // the program has been called with not enough arguments
  if (nArgs <= 2) {
    cerr << "A problem must be specified with the -f option.\n";
    showHelpMessage(programName);
    return 1;
  }

  // extract the filename from the arguments
  if (strcmp(vArgs[1], "-f") == 0) {

    // store the filename
    _smpsFile = string(vArgs[2]);
  }
  else {
    cerr << "A problem must be specified with the -f option.\n";
    showHelpMessage(programName);
    return 1;
  }

  // go through the command line arguments
  for (int arg = 3; arg < nArgs && vArgs[arg] != NULL; ++arg) {

    bool found = false;
    for (int i = 0; i < (int) optionList.size(); ++i) {
      if (strcmp(vArgs[arg], optionList[i].name) == 0) {

	found = true;

	// simple option
	if (!optionList[i].argument) {
	  *optionList[i].variable = 1;
	  break;
	}

	// the option requires an argument
	if (arg + 1 < nArgs) {
	  *optionList[i].variable = atoi(vArgs[++arg]);
	  break;
	}
	else {
	  cout << "The option '" << optionList[i].name
	       << "' requires an argument.\n";
	  return 1;
	}
      }
    }

    // the program has been called with an unrecognized option
    if (!found) {
      cerr << "Option '" << vArgs[arg] << "' not recognised.\n";
      showHelpMessage(programName);
      return 1;
    }
  }

  return 0;
}

/** Add a command line option */
void Options::addOption(const char name[], const char usage[],
			int *variable, const bool arg) {

  OptionAtom newopt = {name, usage, variable, arg};
  optionList.push_back(newopt);
}

/** Print a usage message */
void Options::showHelpMessage(const char* programName) {

  printf("Usage: %s -f <problem.smps> [options]\n"
	 "   where <problem.smps> is the input SMPS file.\n"
	 "Available options:\n", programName);
  for (int i = 0; i < (int) optionList.size(); ++i) {
    printf(" %s   %s\n", optionList[i].name, optionList[i].usage);
  }
}
