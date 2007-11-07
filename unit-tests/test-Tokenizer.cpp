/*
 *  test-Tokenizer.cpp
 *
 *  Unit tests for the Tokenizer class.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include "unit-tests.h"
#include "Tokenizer.h"


static int
testTokenizer(string line, string expToken,
	      const int expLength, const bool expMore = false) {

  int rv;

  char *rline = const_cast<char *>(line.c_str());
  Tokenizer tokenLine(rline);

  string token = tokenLine.getStartNextToken();

  rv = checkEqual(token, expToken, "token");
  rv = checkEqual(tokenLine.getLength(), expLength, "length");
  rv = checkEqual(tokenLine.hasMoreTokens(), expMore, "more");

  return rv;
}

static int
testCountTokens(string line, const int expTokens) {

  int rv;

  char *rline = const_cast<char *>(line.c_str());
  Tokenizer tokenLine(rline);

  int nTokens = tokenLine.countTokens();

  rv = checkEqual(nTokens, expTokens, "nTokens");

  return rv;
}

/**
 *  Call the test function for different inputs.
 */
int unitTokenizer(void) {

  testTokenizer("line", "line", 4);
  testTokenizer(" line", "line", 4);
  testTokenizer("  line", "line", 4);
  testTokenizer("  line  ", "line  ", 4);
  testTokenizer("  line\t", "line\t", 4);
  testTokenizer("\tline  ", "line  ", 4);

  testTokenizer("line cont", "line cont", 4, true);
  testTokenizer(" line cont", "line cont", 4, true);
  testTokenizer("  line cont", "line cont", 4, true);
  testTokenizer("  line  cont", "line  cont", 4, true);
  testTokenizer("  line\tcont", "line\tcont", 4, true);
  testTokenizer("\tline  cont", "line  cont", 4, true);

  testCountTokens("one", 1);
  testCountTokens("", 0);
  testCountTokens("two words", 2);
  testCountTokens("three words line", 3);
  testCountTokens("a line\twith tabs", 4);
  testCountTokens(" short line \t", 2);

  return 0;
}
