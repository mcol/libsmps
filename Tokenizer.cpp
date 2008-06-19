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
  line(rLine),
  pos(rLine),
  length(-1) {
}

/** Destructor */
Tokenizer::~Tokenizer() {}

/** Find the starting position of the next token from the line */
char* Tokenizer::getStartNextToken() {

  // find the next non-space character
  while (*pos == ' ' || *pos == '\t')
    ++pos;

  char *oldpos = pos;

  // check if we have reached the end of the line (marked by '\0')
  if (*pos == 0) {
    length = 0;
    return NULL;
  }

  // pos now points to the next non-space character

  // find the position of the first occurrence of a space or a tab
  // in the string pointed by pos
  char *p2 = strpbrk(pos, " \t");

  // the space character was not found
  if (p2 == NULL) {
    length = strlen(pos);
    pos = &(pos[length]);
  }

  // the space character was found
  else {
    length = p2 - pos;
    pos = p2;
    while (*pos == ' ' || *pos == '\t')
      ++pos;
  }

  if (oldpos[length - 1] == '\n')
    --length;

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
