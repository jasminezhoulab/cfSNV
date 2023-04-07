//
//  read_parser.cpp
//  xfread
//
//  Created by Shuo Li on 6/3/21.
//
#include <Rcpp.h>
using namespace Rcpp;

#include "read_parser.hpp"

#include <math.h>
#include <iostream>

#include "pileup_operators.hpp"


using namespace std;


Parser::Parser() {
    max_base_cnt = DEFAULT_BASE_LEN;
    alloc_chars(DEFAULT_BASE_LEN);
    register_operator(m);

    pkg.in_map['+' - MIN_ASCII] = 1;
    pkg.in_map['-' - MIN_ASCII] = 0;
    pkg.del_map['+' - MIN_ASCII] = 0;
    pkg.del_map['-' - MIN_ASCII] = 1;

    pkg.char_map['.' - MIN_ASCII] = 'R';
    pkg.char_map[',' - MIN_ASCII] = 'r';
    pkg.atcg_map['.' - MIN_ASCII] = 0;
    pkg.atcg_map[',' - MIN_ASCII] = 0;

    for (auto key : "ACGTacgt") {
        pkg.char_map[key - MIN_ASCII] = key;
        pkg.atcg_map[key - MIN_ASCII] = 1;
    }
}


Parser::~Parser() {
    clean_chars();
}


void Parser::register_operator(OpFunction m[MAX_ASCII - MIN_ASCII + 1]) {
    // Set default mapping
    for (auto i = MIN_ASCII; i <= MAX_ASCII; ++i) {
        m['.' - MIN_ASCII] = &match_on_undefined;
    }

    // Define valid operations
    m['.' - MIN_ASCII] = &match_on_forward_reverse;
    m[',' - MIN_ASCII] = &match_on_forward_reverse;
    for (auto key : "ACGTacgt") {  // NOTE: There's an explicit '\0' in this string
        m[key - MIN_ASCII] = &match_on_forward_reverse;
    }
    m['N' - MIN_ASCII] = &mismatch_on_forward_reverse;
    m['n' - MIN_ASCII] = &mismatch_on_forward_reverse;

    m['*' - MIN_ASCII] = &match_on_star;
    m['^' - MIN_ASCII] = &match_on_caret;
    m['$' - MIN_ASCII] = &match_on_dollar;

    m['+' - MIN_ASCII] = &match_on_indel;
    m['-' - MIN_ASCII] = &match_on_indel;
}


void Parser::set_input(int cnt, string &str_base, string &str_qual, string &str_map) {
    if (cnt > max_base_cnt) {
        clean_chars();
        alloc_chars(cnt);
        max_base_cnt = cnt;
    }

    pkg.base_cnt = cnt;

    pkg.base = str_base;
    pkg.qual = str_qual;
    pkg.map = str_map;

    pkg.num_add = 0;
    pkg.num_remove = 0;
    pkg.base_idx = 0;
    pkg.qual_map_idx = 0;
    pkg.ret_idx = 0;
    pkg.ret_pass_idx = 0;

    pkg.num_variant = 0;
    pkg.num_variant_pass = 0;

    for (auto c : "ATCGatcgRr") {  // NOTE: There's an explicit '\0' in this string
        pkg.nuc_cnt_map[c - MIN_ASCII] = 0;
        pkg.nuc_pass_cnt_map[c - MIN_ASCII] = 0;
    }
    return;
}


void Parser::alloc_chars(int len) {
    pkg.ret_mod_base = new char[len+1]();
    pkg.ret_mod_qual = new char[len+1]();
    pkg.ret_mod_map = new char[len+1]();
    pkg.ret_mod_base_pass = new char[len+1]();
    pkg.ret_mod_qual_pass = new char[len+1]();
    pkg.ret_mod_map_pass = new char[len+1]();
}


void Parser::clean_chars() {
    delete pkg.ret_mod_base;
    delete pkg.ret_mod_qual;
    delete pkg.ret_mod_map;
    delete pkg.ret_mod_base_pass;
    delete pkg.ret_mod_qual_pass;
    delete pkg.ret_mod_map_pass;
}


void Parser::parse(Rebase &r) {
    if (pkg.base_cnt != 0) {
        while (pkg.base_idx < pkg.base.length()) {
            // call operator, op = pkg.base[pkg.base_idx]
            m[pkg.base[pkg.base_idx] - MIN_ASCII](pkg);
        }
    }

#if DEBUG_MODE
    translate_to_py1_convention(r);
#else
    r.str_ret_mod_base = string(pkg.ret_mod_base, pkg.ret_idx);
    r.str_ret_mod_qual = string(pkg.ret_mod_qual, pkg.ret_idx);
    r.str_ret_mod_map = string(pkg.ret_mod_map, pkg.ret_idx);
    r.str_ret_mod_base_pass = string(pkg.ret_mod_base_pass, pkg.ret_pass_idx);
    r.str_ret_mod_qual_pass = string(pkg.ret_mod_qual_pass, pkg.ret_pass_idx);
    r.str_ret_mod_map_pass = string(pkg.ret_mod_map_pass, pkg.ret_pass_idx);
#endif

    r.num_add = pkg.num_add;
    r.num_remove = pkg.num_remove;
    r.num_variant = pkg.num_variant;
    r.num_variant_pass = pkg.num_variant_pass;
    r.num_count = (int)r.str_ret_mod_base.length() + r.num_add + r.num_remove;

    for (auto c : "ATCGatcgRr") {  // NOTE: There's an explicit '\0' in this string
        r.nuc_cnt_map[c] = pkg.nuc_cnt_map[c - MIN_ASCII];
        r.nuc_pass_cnt_map[c] = pkg.nuc_pass_cnt_map[c - MIN_ASCII];
    }
}


void Parser::reorder_according_to_base(string &base, string &qual, string &map, string &ret_base, string &ret_qual, string &ret_map) {
    ret_base = "";
    ret_qual = "";
    ret_map = "";

    string temp_a_b = "";
    string temp_c_b = "";
    string temp_g_b = "";
    string temp_t_b = "";

    string temp_a_q = "";
    string temp_c_q = "";
    string temp_g_q = "";
    string temp_t_q = "";

    string temp_a_m = "";
    string temp_c_m = "";
    string temp_g_m = "";
    string temp_t_m = "";
    int idx = 0;
    for (auto c : base) {
        switch (c) {
            case 'R':
            case 'r':
                ret_base += c;
                ret_qual += qual[idx];
                ret_map += map[idx];
                break;
            case 'A':
            case 'a':
                temp_a_b += c;
                temp_a_q += qual[idx];
                temp_a_m += map[idx];
                break;
            case 'C':
            case 'c':
                temp_c_b += c;
                temp_c_q += qual[idx];
                temp_c_m += map[idx];
                break;
            case 'G':
            case 'g':
                temp_g_b += c;
                temp_g_q += qual[idx];
                temp_g_m += map[idx];
                break;
            case 'T':
            case 't':
                temp_t_b += c;
                temp_t_q += qual[idx];
                temp_t_m += map[idx];
                break;
            default:
                cout << "ERROR: Unexpect char " << c << " in base string" << endl;
                exit(5);
                break;
        }
        ++idx;
    }

    ret_base = ret_base + temp_a_b + temp_c_b + temp_g_b + temp_t_b;
    ret_qual = ret_qual + temp_a_q + temp_c_q + temp_g_q + temp_t_q;
    ret_map = ret_map + temp_a_m + temp_c_m + temp_g_m + temp_t_m;
}


void Parser::translate_to_py1_convention(Rebase &r) {
    string temp_str_ret_mod_base = string(pkg.ret_mod_base, pkg.ret_idx);
    string temp_str_ret_mod_qual = string(pkg.ret_mod_qual, pkg.ret_idx);
    string temp_str_ret_mod_map = string(pkg.ret_mod_map, pkg.ret_idx);

    string temp_str_ret_mod_base_pass = string(pkg.ret_mod_base_pass, pkg.ret_pass_idx);
    string temp_str_ret_mod_qual_pass = string(pkg.ret_mod_qual_pass, pkg.ret_pass_idx);
    string temp_str_ret_mod_map_pass = string(pkg.ret_mod_map_pass, pkg.ret_pass_idx);

    // Reorder base string according to R/r, A/a, C/c, G/g, T/t
    reorder_according_to_base(temp_str_ret_mod_base, temp_str_ret_mod_qual, temp_str_ret_mod_map, r.str_ret_mod_base, r.str_ret_mod_qual, r.str_ret_mod_map);
    reorder_according_to_base(temp_str_ret_mod_base_pass, temp_str_ret_mod_qual_pass, temp_str_ret_mod_map_pass, r.str_ret_mod_base_pass, r.str_ret_mod_qual_pass, r.str_ret_mod_map_pass);
}
