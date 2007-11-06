/*
 *  Tokenizer.h
 *
 *  Facilities to tokenize a line of text.
 *
 *  Marco Colombo
 *  School of Mathematics
 *  University of Edinburgh
 *
 */

#ifndef _TOKENIZER_H_
#define _TOKENIZER_H_


/** The Tokenizer class */
class Tokenizer {

 public:

  /** Constructor */
  Tokenizer(char *line);

  /** Destructor */
  ~Tokenizer();

  /** Retrieve the length of the next token */
  int getLength() const {
    return length;
  }

  /** Find the starting position of the next token from the line */
  char* getStartNextToken();

  /** Check if there are other non-space characters in the line */
  bool hasMoreTokens();

 private:

  /** The line to tokenize */
  char *line;

  /** Current position in the line (part of the line still to tokenize) */
  char *pos;

  /** Length of the token extracted */
  int length;

};


#endif /* _TOKENIZER_H_ */