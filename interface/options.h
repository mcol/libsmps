/**
 *  @file options.h
 *
 *  Command line options.
 *
 *  @author Marco Colombo,
 *  School of Mathematics,
 *  University of Edinburgh.
 *
 */

#ifndef _SMPS_OPTIONS_H_
#define _SMPS_OPTIONS_H_

#include <string>
#include <vector>
using namespace std;


/** Command line options for a generic SMPS interface */
class Options {

  /** Number of arguments in the vArgs array */
  const int nArgs;

  /** Values of the command line arguments */
  const char **vArgs;

  /** Whether the problem should not be solved */
  int _dontSolve;

  /** Whether the output has to be redirected to file */
  int _outputToFile;

  /** Whether an mps file has to be produced  */
  int _writeMps;

  /** Whether the solution has to be printed to screen */
  int _printSolution;

  /** Name of the smps file */
  string smpsFile;

  /** Definition of each option */
  struct OptionAtom {

    /** The way of calling the option */
    const char *name;

    /** Description of the option for the usage string */
    const char *usage;

    /** Pointer to the variable to set */
    int *variable;

    /** Wheter the option requires an argument */
    bool argument;

  };

  /** The list of options to be parsed */
  vector<OptionAtom> optionList;

  /** Print an usage message */
  void showHelpMessage(const char *programName);

 protected:

  /** Add a command line option */
  void addOption(const char name[], const char description[],
		 int *variable, const bool arg = false);

 public:

  /** Constructor */
  Options(const int argc = 0, const char *argv[] = NULL);

  /** Destructor */
  ~Options();

  /** Parse the command line options */
  int parse(void);

  /** Retrieve the name of the smps file */
  string getSmpsFile(void) const {
    return smpsFile;
  }

  /** Retrieve the value of the dontSolve option */
  int dontSolve(void) const {
    return _dontSolve;
  }

  /** Retrieve the value of the outputToFile option */
  int outputToFile(void) const {
    return _outputToFile;
  }

  /** Retrieve the value of the writeMps option */
  int writeMps(void) const {
    return _writeMps;
  }

  /** Retrieve the value of the printSolution option */
  int printSolution(void) const {
    return _printSolution;
  }

};

#endif /* _SMPS_OPTIONS_H_ */
