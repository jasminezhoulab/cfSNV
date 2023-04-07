//
//  read_parser.hpp
//  xfread
//
//  Created by Shuo Li on 6/3/21.
//

#ifndef read_parser_hpp
#define read_parser_hpp

#include <Rcpp.h>
using namespace Rcpp;

#include <stdio.h>
#include <unordered_map>
#include <map>
#include <string>

#define DEBUG_MODE 1

#define DEFAULT_BASE_LEN 1000
#define DEFAULT_BASE_MULT 10
#define BASEQUAL_THRESHOLD 0.05
#define MAPQUAL_THRESHOLD 0.3
#define COUNT_VAR_TUMOR_HIGHQUAL 3
#define TRIALLELE_VAF_RATIO 0.5
#define TRIALLELE_VAF 0.02
#define DEPTH_FOR_DETECTION_TUMOR 4
#define DEPTH_FOR_DETECTION_NORMAL 4
#define NORMAL_COUNT_ALT 7

#define SOMATIC_VAF_THRESHOLD_IN_NORMAL 0.12
#define GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER 2

#define MAX_ASCII 'z'
#define MIN_ASCII '\0'
#define READ_BUF_SIZE 256*1024


typedef struct {
    std::string str_ret_mod_base;  // one on one translation of input base line (filter out special chars)
    std::string str_ret_mod_qual;  // the same as the input qual line
    std::string str_ret_mod_map;  // the same as the input map line

    std::string str_ret_mod_base_pass;  // filter based on qual AND map conditions
    std::string str_ret_mod_qual_pass;  // corresponding qual string after base_pass filtering
    std::string str_ret_mod_map_pass;  // corresponding map string after base_pass filtering

    // Additional data for validating if the input line should be included in the output file
    int num_add;
    int num_remove;
    int num_variant;
    int num_variant_pass;

    int num_count;

    std::unordered_map<char, int> nuc_cnt_map;
    std::unordered_map<char, int> nuc_pass_cnt_map;
} Rebase;


typedef struct {
    // constants
    int in_map['.' - MIN_ASCII + 1];
    int del_map['.' - MIN_ASCII + 1];
    int atcg_map[MAX_ASCII - MIN_ASCII + 1];
    char char_map[MAX_ASCII - MIN_ASCII + 1];

    // buffers
    char *ret_mod_base;  // len <= base_cnt
    char *ret_mod_qual;  // len <= base_cnt
    char *ret_mod_map;  // len <= base_cnt
    char *ret_mod_base_pass;  // len <= base_cnt
    char *ret_mod_qual_pass;  // len <= base_cnt
    char *ret_mod_map_pass;  // len <= base_cnt

    // inputs
    int base_cnt;

    std::string base;
    std::string qual;
    std::string map;

    // temp variables
    int base_idx;
    int qual_map_idx;
    int ret_idx;
    int ret_pass_idx;

    // results
    int num_add;
    int num_remove;

    int num_variant;
    int num_variant_pass;

    int nuc_cnt_map[MAX_ASCII - MIN_ASCII + 1];
    int nuc_pass_cnt_map[MAX_ASCII - MIN_ASCII + 1];
} ReadPkg;


typedef void (*OpFunction)(ReadPkg &);  // function pointer type


class Parser{
public:
    Parser();
    ~Parser();

    void set_input(int cnt, std::string &str_base, std::string &str_qual, std::string &str_map);
    void parse(Rebase &r);


private:
    int max_base_cnt;
    OpFunction m[MAX_ASCII - MIN_ASCII + 1];

    // Temporary buffer
    ReadPkg pkg;

    std::string str_ret_mod_base;
    std::string str_ret_mod_qual;
    std::string str_ret_mod_map;
    std::string str_ret_mod_base_pass;
    std::string str_ret_mod_qual_pass;
    std::string str_ret_mod_map_pass;

    void alloc_chars(int len);
    void clean_chars();

    void register_operator(OpFunction m[MAX_ASCII - MIN_ASCII + 1]);
    void translate_to_py1_convention(Rebase &r);
    void reorder_according_to_base(std::string &base, std::string &qual, std::string &map, std::string &ret_base, std::string &ret_qual, std::string &ret_map);
};


#endif /* read_parser_hpp */
