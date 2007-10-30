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

#include <stdio.h>
#include <cstring>
#include "unit-tests.h"
#include "Tokenizer.h"

extern int nTests;
extern int nFails;

int testTokenizer(char *line, char *expToken,
		  const int expLength, const bool expMore = false) {

  int rv = TEST_SUCCESS;

  Tokenizer tokenLine(line);

  char *token = tokenLine.getStartNextToken();
  nTests++;

  int length = tokenLine.getLength();
  nTests++;

  bool more = tokenLine.hasMoreTokens();
  nTests++;

  if (strcmp(token, expToken) != 0) {
    printf(" | FAIL: token: '%s' (exp: '%s')\n", token, expToken);
    nFails++;
    rv = TEST_FAILURE;
  }

  if (length != expLength) {
    printf(" * Testing line: '%s'\n", line);
    printf(" | FAIL: length: %d (exp: %d)\n", length, expLength);
    nFails++;
    rv = TEST_FAILURE;
  }

  if (more != expMore) {
    printf(" * Testing line: '%s'\n", line);
    printf(" | FAIL: more: %d (exp: %d)\n", more, expMore);
    nFails++;
    rv = TEST_FAILURE;
  }

  return 0;
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
