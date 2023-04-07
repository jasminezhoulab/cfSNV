//
//  pileup_parser.cpp
//  xfread
//
//  Created by Shuo Li on 6/5/21.
//
#include <Rcpp.h>
using namespace Rcpp;

#include "pileup_parser.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/math/distributions/binomial.hpp>

#include "read_parser.hpp"


using namespace boost::math;
using namespace std;


std::vector<std::string>
split(std::string const& original, char separator) {
    std::vector<std::string> results;
    std::string::const_iterator start = original.begin();
    std::string::const_iterator end = original.end();
    std::string::const_iterator next = std::find(start, end, separator);
    while (next != end) {
        results.push_back(std::string(start, next));
        start = next + 1;
        next = std::find(start, end, separator);
    }
    results.push_back(std::string(start, next));
    return results;
}


char find_major_variant(unordered_map<char, int> &nc_cnt_map, unordered_map<char, int> &ef_cnt_map) {
    int arr_cnt[4] = {0, 0, 0, 0};
    char arr_atcg[4] = {'A', 'C', 'T', 'G'};

    arr_cnt[0] = nc_cnt_map['a'] + nc_cnt_map['A'] + 2*ef_cnt_map['a'] + 2*ef_cnt_map['A'];
    arr_cnt[1] = nc_cnt_map['c'] + nc_cnt_map['C'] + 2*ef_cnt_map['c'] + 2*ef_cnt_map['C'];
    arr_cnt[2] = nc_cnt_map['t'] + nc_cnt_map['T'] + 2*ef_cnt_map['t'] + 2*ef_cnt_map['T'];
    arr_cnt[3] = nc_cnt_map['g'] + nc_cnt_map['G'] + 2*ef_cnt_map['g'] + 2*ef_cnt_map['G'];

    int max_variant_cnt = 0;
    char major_variant = '\0';
    for (int i=3; i>=0; --i) {
        if (arr_cnt[i] >= max_variant_cnt) {
            major_variant = arr_atcg[i];
            max_variant_cnt = arr_cnt[i];
        }
    }

    return major_variant;
}


int get_germline_variant_count_threshold(float depth) {
    return int(log10(depth))-1;
}


float get_normal_count_binom_threshold(float depth) {
    if (depth < 1000.0)
        return 0.05;
    else
        return 0.2;
}


bool filter_triallelic_position(unordered_map<char, int> &cnt_map, int base_len, char major_variant, float depth) {
    int arr_cnt[4] = {0, 0, 0, 0};
    char arr_atcg[4] = {'A', 'T', 'C', 'G'};

    arr_cnt[0] = cnt_map['a'] + cnt_map['A'];
    arr_cnt[1] = cnt_map['t'] + cnt_map['T'];
    arr_cnt[2] = cnt_map['c'] + cnt_map['C'];
    arr_cnt[3] = cnt_map['g'] + cnt_map['G'];

    int max_other_variant_cnt = 0;
    char other_variant = '\0';
    for (int i=0; i<4; ++i) {
        if (arr_atcg[i] == major_variant) continue;
        if (arr_cnt[i] >= max_other_variant_cnt) {
            other_variant = arr_atcg[i];
            max_other_variant_cnt = arr_cnt[i];
        }
    }

    float other_variant_f = float(max_other_variant_cnt) / float(base_len);
    int major_variant_cnt = cnt_map[major_variant] + cnt_map[tolower(major_variant)];
    float major_variant_f = float(major_variant_cnt) / float(base_len);

    int triallele_cnt = int(log10(depth)) + 2;


    if (major_variant_cnt == 0 || other_variant_f*1.0/major_variant_f > TRIALLELE_VAF_RATIO)
        return false;
    if (other_variant_f > TRIALLELE_VAF && base_len >= 80)
        return false;
    if (max_other_variant_cnt > triallele_cnt && base_len <= 140)
        return false;
    return true;
}

// [[Rcpp::export]]
void ParseFile(std::string in_file_name, std::string out_file_name, std::string indel_file_name, float depth, int MIN_PASS_SUPPORT_COUNT) {
    const size_t bufsize = READ_BUF_SIZE;
    char buf[bufsize];

    ifstream in_file;
    in_file.rdbuf()->pubsetbuf(buf, bufsize);

    ofstream out_file, indel_file;
    in_file.open(in_file_name);

    Parser p;
    Rebase tumor;
    Rebase normal;
    Rebase extendedFrags;
    Rebase notCombined ;
    int line_cnt = 0;

    if (in_file.is_open()){
        out_file.open(out_file_name);
        if (!out_file.is_open()) {
            cout << "Cannot write to " << out_file_name << endl;
            exit(1);
        }
        indel_file.open(indel_file_name);
        if (!indel_file.is_open()) {
            cout << "Cannot write to " << indel_file_name << endl;
            exit(1);
        }

        string str_temp;
        while(getline(in_file, str_temp)) {
            line_cnt += 1;
            auto sp = split(str_temp, '\t');

            if (sp[3] == "0") {
                continue;
            }
            p.set_input(stoi(sp[3]), sp[4], sp[5], sp[6]);
            p.parse(tumor);

            if (sp[7] == "0") {
                continue;
            }
            p.set_input(stoi(sp[7]), sp[8], sp[9], sp[10]);
            p.parse(normal);


            p.set_input(stoi(sp[11]), sp[12], sp[13], sp[14]);
            p.parse(extendedFrags);
            p.set_input(stoi(sp[15]), sp[16], sp[17], sp[18]);
            p.parse(notCombined);


            // Decide if this input line should be included in output file
            if (normal.num_remove >=2 || tumor.num_remove >= 2 || extendedFrags.num_remove >= 2 || notCombined.num_remove >= 2) {
                indel_file << (sp[0] + '\t' + sp[1] + '\t' + "delete" + '\t' + to_string(normal.num_remove) + '\t' + to_string(normal.num_count) + '\t' + to_string(tumor.num_remove) + '\t' + to_string(tumor.num_count) + '\t' + to_string(extendedFrags.num_remove) + '\t' + to_string(extendedFrags.num_count) + '\t' + to_string(notCombined.num_remove) + '\t' + to_string(notCombined.num_count) + '\n');
                continue;
            }
            if (normal.num_add >= 2 || tumor.num_add >= 2 || extendedFrags.num_add >= 2 || notCombined.num_add >=2) {
                indel_file << (sp[0] + '\t' + sp[1] + '\t' + "insert" + '\t' + to_string(normal.num_add) + '\t' + to_string(normal.num_count) + '\t' + to_string(tumor.num_add) + '\t' + to_string(tumor.num_count) + '\t' + to_string(extendedFrags.num_add) + '\t' + to_string(extendedFrags.num_count) + '\t' + to_string(notCombined.num_add) + '\t' + to_string(notCombined.num_count) + '\n');
                continue;
            }


            if (tumor.num_variant_pass == 0
                && (extendedFrags.num_variant_pass + notCombined.num_variant_pass == 0)) {
                continue;
            }
            if (tumor.str_ret_mod_base_pass.length() == 0
                || normal.str_ret_mod_base_pass.length() == 0
                ) {
                continue;
            }

            auto major_variant_base = find_major_variant(notCombined.nuc_cnt_map, extendedFrags.nuc_cnt_map);
            auto major_variant_base_pass = find_major_variant(notCombined.nuc_pass_cnt_map, extendedFrags.nuc_pass_cnt_map);

            if (major_variant_base != major_variant_base_pass){
                continue;
            }


            auto triallele_tumor = filter_triallelic_position(tumor.nuc_cnt_map, (int)tumor.str_ret_mod_base.length(), major_variant_base, depth);
            auto triallele_tumor_pass = filter_triallelic_position(tumor.nuc_pass_cnt_map, (int)tumor.str_ret_mod_base_pass.length(), major_variant_base_pass, depth);

            if ((! triallele_tumor) || (! triallele_tumor_pass)){
                continue;
            }

            auto cnt_var_normal_pass = normal.nuc_pass_cnt_map[major_variant_base_pass] + normal.nuc_pass_cnt_map[tolower(major_variant_base_pass)];

            auto cnt_alt_normal_pass = normal.str_ret_mod_base_pass.length() - normal.nuc_pass_cnt_map['r'] - normal.nuc_pass_cnt_map['R'];

            auto cnt_normal_pass = normal.str_ret_mod_base_pass.length();

            auto cnt_tumor_merge_pass = extendedFrags.str_ret_mod_base_pass.length() + notCombined.str_ret_mod_base_pass.length();
            if ((cnt_tumor_merge_pass < DEPTH_FOR_DETECTION_TUMOR)
                ||
                (tumor.str_ret_mod_base_pass.length() < DEPTH_FOR_DETECTION_TUMOR)
                ){
                continue;
            }

            if (cnt_normal_pass < DEPTH_FOR_DETECTION_NORMAL){
                continue;
            }

            int NORMAL_COUNT_VAR = floor(0.5*MIN_PASS_SUPPORT_COUNT);
            if (cnt_var_normal_pass > NORMAL_COUNT_VAR
                || cnt_alt_normal_pass > NORMAL_COUNT_ALT){
                continue;
            }

            auto cnt_var_tumor_pass = tumor.nuc_pass_cnt_map[major_variant_base_pass] + tumor.nuc_pass_cnt_map[tolower(major_variant_base_pass)];

            float vaf_normal_pass = float(cnt_var_normal_pass)/float(cnt_normal_pass);

            auto cnt_var_tumor_merge_pass = notCombined.nuc_pass_cnt_map[major_variant_base_pass] + notCombined.nuc_pass_cnt_map[tolower(major_variant_base_pass)] + extendedFrags.nuc_pass_cnt_map[major_variant_base_pass] + extendedFrags.nuc_pass_cnt_map[tolower(major_variant_base_pass)];

            float vaf_tumor_merge_pass = float(cnt_var_tumor_merge_pass)/float(cnt_tumor_merge_pass);

            float vaf_tumor_pass = float(cnt_var_tumor_pass)/float(tumor.str_ret_mod_base_pass.length());

            if (vaf_normal_pass > SOMATIC_VAF_THRESHOLD_IN_NORMAL){
                continue;
            }

            auto cnt_var_normal = normal.nuc_cnt_map[major_variant_base_pass] + normal.nuc_cnt_map[tolower(major_variant_base_pass)];

            if (cnt_var_normal > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER
                || cnt_var_normal_pass > get_germline_variant_count_threshold(depth)){
                continue;
            }

            if (cnt_var_normal_pass > 0
                && (vaf_normal_pass >= vaf_tumor_pass
                    || vaf_normal_pass >= vaf_tumor_merge_pass
                    || cnt_var_normal_pass >= cnt_var_tumor_pass
                    || cnt_var_normal_pass >= cnt_var_tumor_merge_pass)){
                continue;
            }

            if(cnt_var_normal_pass > 0 &&
               cdf(binomial((int)normal.str_ret_mod_base_pass.length(), vaf_tumor_pass), cnt_var_normal_pass)
                > get_normal_count_binom_threshold(depth)
               ){
                continue;
            }


            auto cnt_alt_tumor_pass = tumor.str_ret_mod_base_pass.length() - tumor.nuc_pass_cnt_map['r'] - tumor.nuc_pass_cnt_map['R'];
            auto cnt_alt_tumor_merge_pass = notCombined.str_ret_mod_base_pass.length() - notCombined.nuc_pass_cnt_map['r'] - notCombined.nuc_pass_cnt_map['R'] + extendedFrags.str_ret_mod_base_pass.length() - extendedFrags.nuc_pass_cnt_map['r'] - extendedFrags.nuc_pass_cnt_map['R'];

            if (cnt_alt_tumor_pass < COUNT_VAR_TUMOR_HIGHQUAL
                && cnt_alt_tumor_merge_pass < COUNT_VAR_TUMOR_HIGHQUAL){
                continue;
            }

            out_file << (sp[0] + '\t' + sp[1] + '\t' + (char)toupper(sp[2][0]) + '\t' + sp[3] + '\t' + tumor.str_ret_mod_base + '\t' + tumor.str_ret_mod_qual + '\t' + tumor.str_ret_mod_map + '\t' + sp[7] + '\t' + normal.str_ret_mod_base + '\t' + normal.str_ret_mod_qual + '\t' + normal.str_ret_mod_map + '\t' + sp[11] + '\t' + extendedFrags.str_ret_mod_base + '\t' + extendedFrags.str_ret_mod_qual + '\t' + extendedFrags.str_ret_mod_map + '\t' + sp[15] + '\t' + notCombined.str_ret_mod_base + '\t' + notCombined.str_ret_mod_qual + '\t' + notCombined.str_ret_mod_map + '\t' + tumor.str_ret_mod_base_pass + '\t' + tumor.str_ret_mod_qual_pass + '\t' + tumor.str_ret_mod_map_pass + '\t' + normal.str_ret_mod_base_pass + '\t' + normal.str_ret_mod_qual_pass + '\t' + normal.str_ret_mod_map_pass + '\t' + extendedFrags.str_ret_mod_base_pass + '\t' + extendedFrags.str_ret_mod_qual_pass + '\t' + extendedFrags.str_ret_mod_map_pass + '\t' + notCombined.str_ret_mod_base_pass + '\t' + notCombined.str_ret_mod_qual_pass + '\t' + notCombined.str_ret_mod_map_pass + '\n');
        }

        in_file.close();
        out_file.close();
        indel_file.close();
    } else {
        cout << "Cannot read " << in_file_name << endl;
        exit(1);
    }
}
