/*
 *  Tokenizer.cpp
 *
 *  Facilities to tokenize a line of text.
 *
 *  Andreas Grothey
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#include <stdlib.h>
#include <cstring>
#include "Tokenizer.h"


/** Constructor */
Tokenizer::Tokenizer(char *rLine) :
  pos(rLine),
  ready(false),
  length(-1) {
}

/** Destructor */
Tokenizer::~Tokenizer() {}

/** Extract a token from the line */
char* Tokenizer::getToken() {

  // token delimiters
  const char* blanks = " \t";

  char *param = NULL;

  // strtok needs to be called with the line to tokenize the first time,
  // and with NULL all subsequent times
  if (!ready) {

    // start tokenizing from the current position
    strcpy(line, pos);
    param = line;
    ready = true;
  }

  // extract the token
  char *token = strtok(param, blanks);

  // compute the length of the token
  length = token ? strlen(token) : 0;

  // find the next non-space character
  while (*pos == ' ' || *pos == '\t')
    ++pos;

  // move the pointer to the end of the token
  pos += length;

  return token;
}

/** Find the starting position of the next token from the line */
char* Tokenizer::getStartNextToken() {

  // find the next non-space character
  while (*pos == ' ' || *pos == '\t')
    ++pos;

  // since pos has changed, we will have to reset the tokenizer
  ready = false;

  // check if we have reached the end of the line (marked by '\0')
  if (*pos == 0) {
    return NULL;
  }

  // pos points to the next non-space character (beginning of the next token)
  char *oldpos = pos;

  // now we move pos to the end of the next token

  // find the position of the first occurrence of a space or a tab
  // in the string pointed by pos
  char *p2 = strpbrk(pos, " \t");

  // the space character was not found
  if (p2 == NULL) {
    int len = strlen(pos);
    pos = &(pos[len]);
  }

  // the space character was found
  else {
    pos = p2;
    while (*pos == ' ' || *pos == '\t')
      ++pos;
  }

  return oldpos;
}

/** Check if there are other non-space characters in the line */
bool Tokenizer::hasMoreTokens() {

  if (!pos)
    return false;

  // skip all whitespace characters
  while (*pos == ' ' || *pos == '\t' || *pos == '\n')
    ++pos;

  // the end of the line has been reached
  if (*pos == 0)
    return false;

  return true;
}

/** Count the number of tokens in the string */
int Tokenizer::countTokens() {

  int nTokens = 0;

  // store the current position
  char *oldpos = pos;

  while (hasMoreTokens()) {
    getStartNextToken();
    ++nTokens;
  }

  // reset the position
  pos = oldpos;

  return nTokens;
}
