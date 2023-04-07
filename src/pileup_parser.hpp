//
//  pileup_parser.hpp
//  xfread
//
//  Created by Shuo Li on 6/5/21.
//

#ifndef pileup_parser_hpp
#define pileup_parser_hpp

#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>
#include <string>

void ParseFile(std::string in_file_name, std::string out_file_name, std::string indel_file_name, float depth);

#endif /* pileup_parser_hpp */
