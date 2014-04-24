///////////////////////////////////////////////////////
//
//  Parser.h
//
//  Purpose:  Conatins generalized parsing tools, such as tokenizer
//            parsers, and lexer
//       
///////////////////////////////////////////////////////


#ifndef PARSER_H_
#define PARSER_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>



namespace Utilities
{
  
  using std::cerr;
  using std::endl;
  using std::vector;
  using std::string;
  using std::ifstream;
  typedef size_t Size_Type;  
  
  void Tokenize( vector<vector<string> > &vsTokens, const string &sBuf,
                 const string &sDelimiters);

}

#endif
