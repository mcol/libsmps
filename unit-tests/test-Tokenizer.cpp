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
testTokenizer(const string line, const string expToken,
	      const int expLength, const bool expMore = false) {

  int rv;

  char *rline = const_cast<char *>(line.c_str());
  Tokenizer tokenLine(rline);

  char *rtoken = tokenLine.getToken();
  string token = rtoken ? rtoken : "";

  rv = checkEqual(token, expToken, "token");
  rv = checkEqual(tokenLine.getLength(), expLength, "length");
  rv = checkEqual(tokenLine.hasMoreTokens(), expMore, "more");

  return rv;
}

static int
testGetStartNextToken(const string line, const string expLine) {

  int rv;

  char *rline = const_cast<char *>(line.c_str());
  Tokenizer tokenLine(rline);

  char *rtoken = tokenLine.getStartNextToken();
  string token = rtoken ? rtoken : "";

  rv = checkEqual(token, expLine, "nextToken");

  return rv;
}

static int
testMultipleTokens(const string line, const string expToken,
		   const string expToken2 = " ") {

  int rv;

  char *rline = const_cast<char *>(line.c_str());
  Tokenizer tokenLine(rline);

  tokenLine.getStartNextToken();
  char *rtoken = tokenLine.getToken();
  string token = rtoken ? rtoken : "";

  rv = checkEqual(token, expToken, "token");

  if (expToken2 != " ") {
    rtoken = tokenLine.getToken();
    token  = rtoken ? rtoken : "";
    rv = checkEqual(token, expToken2, "token2");
  }

  return rv;
}

static int
testCountTokens(const string line, const int expTokens) {

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

  setFamily("testTokenizer");
  testTokenizer("", "", 0);
  testTokenizer("line", "line", 4);
  testTokenizer(" line", "line", 4);
  testTokenizer("  line", "line", 4);
  testTokenizer("  line  ", "line", 4);
  testTokenizer("  line\t", "line", 4);
  testTokenizer("\tline  ", "line", 4);
  testTokenizer("line cont", "line", 4, true);
  testTokenizer(" line cont", "line", 4, true);
  testTokenizer("  line cont", "line", 4, true);
  testTokenizer("  line  cont", "line", 4, true);
  testTokenizer("  line\tcont", "line", 4, true);
  testTokenizer("\tline  cont", "line", 4, true);

  setFamily("testGetStartNextToken");
  testGetStartNextToken("", "");
  testGetStartNextToken("line", "line");
  testGetStartNextToken(" line", "line");
  testGetStartNextToken("  line", "line");
  testGetStartNextToken("  line  ", "line  ");
  testGetStartNextToken("  line\t", "line\t");
  testGetStartNextToken("\tline  ", "line  ");
  testGetStartNextToken("line cont", "line cont");
  testGetStartNextToken(" line cont", "line cont");
  testGetStartNextToken("  line cont", "line cont");
  testGetStartNextToken("  line  cont", "line  cont");
  testGetStartNextToken("  line\tcont", "line\tcont");
  testGetStartNextToken("\tline  cont", "line  cont");

  setFamily("testMultipleTokens");
  testMultipleTokens("", "");
  testMultipleTokens("first ", "");
  testMultipleTokens(" first second", "second");
  testMultipleTokens("  first second\t ", "second");
  testMultipleTokens("  first  second third", "second");
  testMultipleTokens("  first\tsecond third", "second");
  testMultipleTokens("\tfirst  second\tthird  ", "second");
  testMultipleTokens("", "", "");
  testMultipleTokens(" first second", "second", "");
  testMultipleTokens("  first  second third", "second", "third");
  testMultipleTokens("  first\tsecond third\t", "second", "third");
  testMultipleTokens("\tfirst  second\tthird  ", "second", "third");
  testMultipleTokens(" first second\t third fourth", "second", "third");

  setFamily("testCountTokens");
  testCountTokens("one", 1);
  testCountTokens("", 0);
  testCountTokens("two words", 2);
  testCountTokens("three words line", 3);
  testCountTokens("a line\twith tabs", 4);
  testCountTokens("another line \t with tabs\n", 4);
  testCountTokens(" short line \t", 2);
  testCountTokens(" short line \0 should not read this", 2);

  return 0;
}
