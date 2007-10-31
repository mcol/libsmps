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

extern int nTests;
extern int nFails;

int testTokenizer(char *line, char *expToken,
		  const int expLength, const bool expMore = false) {

  int rv;

  Tokenizer tokenLine(line);

  char *token = tokenLine.getStartNextToken();

  rv = checkEqual(token, expToken, "token");
  rv = checkEqual(tokenLine.getLength(), expLength, "length");
  rv = checkEqual(tokenLine.hasMoreTokens(), expMore, "more");

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

  return 0;
}
