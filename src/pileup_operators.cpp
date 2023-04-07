//
//  pileup_operators.cpp
//  xfread
//
//  Created by Shuo Li on 6/5/21.
//
#include <Rcpp.h>
using namespace Rcpp;

#include "pileup_operators.hpp"

#include <math.h>
#include <iostream>

using namespace std;

void match_on_forward_reverse(ReadPkg &pkg) {
    auto c = pkg.base[pkg.base_idx] - MIN_ASCII;

    pkg.ret_mod_base[pkg.ret_idx] = pkg.char_map[c];
    pkg.ret_mod_qual[pkg.ret_idx] = pkg.qual[pkg.qual_map_idx];
    pkg.ret_mod_map[pkg.ret_idx] = pkg.map[pkg.qual_map_idx];

    if (
        (pow(10.0, (-float(int(pkg.qual[pkg.qual_map_idx])-33)/10.0)) <= BASEQUAL_THRESHOLD)
        &&
        (pow(10.0, (-float(int(pkg.map[pkg.qual_map_idx])-33)/10.0)) <= MAPQUAL_THRESHOLD)
        ) {

        pkg.ret_mod_base_pass[pkg.ret_pass_idx] = pkg.char_map[c];
        pkg.ret_mod_qual_pass[pkg.ret_pass_idx] = pkg.qual[pkg.qual_map_idx];
        pkg.ret_mod_map_pass[pkg.ret_pass_idx] = pkg.map[pkg.qual_map_idx];
        pkg.ret_pass_idx += 1;

        pkg.num_variant_pass += pkg.atcg_map[c];
        pkg.nuc_pass_cnt_map[int(pkg.char_map[c] - MIN_ASCII)] += 1;
    }

    pkg.qual_map_idx += 1;
    pkg.base_idx += 1;
    pkg.ret_idx += 1;
    pkg.num_variant += pkg.atcg_map[c];
    pkg.nuc_cnt_map[int(pkg.char_map[c] - MIN_ASCII)] += 1;
}


void mismatch_on_forward_reverse(ReadPkg &pkg) {
    pkg.qual_map_idx += 1;
    pkg.base_idx += 1;
}


void match_on_star(ReadPkg &pkg) {
    pkg.qual_map_idx += 1;
    pkg.base_idx += 1;
    pkg.num_remove += 1;
}

void match_on_caret(ReadPkg &pkg) {
    pkg.base_idx += 2;
}


void match_on_dollar(ReadPkg &pkg) {
    pkg.base_idx += 1;
}


void match_on_indel(ReadPkg &pkg) {
    auto op = pkg.base[pkg.base_idx] - MIN_ASCII;
    pkg.num_add += pkg.in_map[op];
    pkg.num_remove += pkg.del_map[op];

    pkg.base_idx += 1;  // skip current +/-

    int skip_num = 0;
    while (
           pkg.base[pkg.base_idx] >= '0' && pkg.base[pkg.base_idx] <= '9'
           ) {

        auto temp_num = (int)(pkg.base[pkg.base_idx] - '0');
        skip_num = skip_num * 10 + temp_num;
        pkg.base_idx += 1;
    }

    pkg.base_idx += skip_num;
}


void match_on_undefined(ReadPkg &pkg) {
    std::cout << "Operator " << pkg.base[pkg.base_idx] << " is NOT defined." << std::endl;
    exit(1);
}

