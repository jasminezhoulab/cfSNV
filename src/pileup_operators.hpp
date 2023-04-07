//
//  pileup_operators.hpp
//  xfread
//
//  Created by Shuo Li on 6/5/21.
//

#ifndef pileup_operators_hpp
#define pileup_operators_hpp

#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>

#include "read_parser.hpp"

void match_on_forward_reverse(ReadPkg &pkg);
void mismatch_on_forward_reverse(ReadPkg &pkg);
void match_on_star(ReadPkg &pkg);
void match_on_caret(ReadPkg &pkg);
void match_on_dollar(ReadPkg &pkg);
void match_on_indel(ReadPkg &pkg);
void match_on_undefined(ReadPkg &pkg);

#endif /* pileup_operators_hpp */
