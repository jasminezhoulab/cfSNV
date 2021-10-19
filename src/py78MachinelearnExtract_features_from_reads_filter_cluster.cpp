//
// Created by Ran Hu on 8/1/21.
//
#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <algorithm>
#include <stdlib.h>
#include <bitset>
#include <cstring>

#include "py78MachinelearnExtract_features_from_reads_filter_cluster.h"
#include "_helper.h"

std::vector<std::string> make_homopolymer(int HOMOPOLYMER_SIZE) {

    char NUCLEOTIDE[4] = {'A', 'C', 'T', 'G'};
    std::vector<std::string> HOMOPOLYMER;
    for (char &l : NUCLEOTIDE)
        HOMOPOLYMER.push_back(std::string(HOMOPOLYMER_SIZE, l));
    return HOMOPOLYMER;
}

std::vector<std::vector<std::string>> read_bed_for_feature(std::string fastabed, int window) {

    std::ifstream file(fastabed);
    std::string lineText;
    std::vector<std::vector<std::string>> preparebed;
    std::string base;

    while (std::getline(file, lineText)) {
        std::vector<std::string> sp = split(lineText, "\t");
        std::vector<std::string> oneLine;
        oneLine.push_back(sp[0]);
        oneLine.push_back(std::to_string(std::stoi(sp[1]) + window));
        oneLine.push_back(std::to_string(std::stoi(sp[2]) - window));
        oneLine.push_back(toUpper(sp[4]));
        base = toUpper(sp[4])[window];
        oneLine.push_back(base);
        oneLine.push_back(sp[3]);
        preparebed.push_back(oneLine);
    }
    return preparebed;
}

void build_mapping_to_var(const std::vector<std::vector<std::string>> &preparebed, std::map<std::string, int> &mapping_to_var,
                          std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<int>>>> &matrix) {

    std::vector<std::vector<std::string>> loc;
    for (auto const &s : preparebed) {
        std::vector<std::string> oneLoc;
        oneLoc.push_back(s[0]);
        oneLoc.push_back(s[2]);
        loc.push_back(oneLoc);
    }
//    int k = 0;
//    for (std::vector<std::string> &s : loc) {
//        for (std::vector<std::string>::iterator i = s.begin(); i != s.end(); ++i) {
//            std::cout << *i << '\t';
//        }
//        std::cout << '\n';
//        k++;
//        if (k == 5)
//            break;
//    }
//    std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;

    std::stable_sort(loc.begin(), loc.end(), cmp1);
//
//    k = 0;
//    for (std::vector<std::string> &s : loc) {
//        for (std::vector<std::string>::iterator i = s.begin(); i != s.end(); ++i) {
//            std::cout << *i << '\t';
//        }
//        std::cout << '\n';
//        k++;
//        if (k == 5)
//            break;
//    }
//    std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;

    std::stable_sort(loc.begin(), loc.end(), cmp0);
//
//    int k = 0;
//    for (std::vector<std::string> &s : loc) {
//        for (std::vector<std::string>::iterator i = s.begin(); i != s.end(); ++i) {
//            std::cout << *i << '\t';
//        }
//        std::cout << '\n';
//        k++;
//        if (k == 100)
//            break;
//    }
//    std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;

    std::vector<std::vector<std::vector<std::string>>> neighbor_list;
    std::vector<std::vector<std::string>> neighbor;
    int n = 0;
    int flag = 0;

    for (int i = 0; i < loc.size()-1; i++) {
        if (loc[i][0] == loc[i+1][0]) {
            if (std::abs(std::stoi(loc[i][1]) - std::stoi(loc[i+1][1])) < 20) {
                n++;
                if (flag == 0) {
                    flag = 1;
                    neighbor.push_back(loc[i]);
                    neighbor.push_back(loc[i+1]);
                } else {
                    neighbor.push_back(loc[i+1]);
                }
            } else {
                if (flag == 1)
                    neighbor_list.push_back(neighbor);
                flag = 0;
                neighbor.clear();
            }
        } else {
            if (flag == 1)
                neighbor_list.push_back(neighbor);
            flag = 0;
            neighbor.clear();
        }
    }
    for (int i = 0; i < neighbor_list.size(); i++) {
        neighbor = neighbor_list[i];
        std::vector<std::string> neighbor_string;
        for (std::vector<std::string> &s : neighbor) {
            std::string one_string = join(std::vector<std::string>(s.begin(), s.end()), "-");
            neighbor_string.push_back(one_string);
            mapping_to_var[one_string] = i;
        }
        int l = neighbor.size();
        std::pair<std::vector<std::string>, std::vector<std::vector<int>>> matrix_element;
        std::vector<std::vector<int>> vec(l, std::vector<int> (l, 0));
        matrix_element.first = neighbor_string;
        matrix_element.second = vec;
        matrix.push_back(matrix_element);
    }
}

//std::string extract_mapping_quality(std::string read) {
//
//    std::string phred = split(read, "\t")[4];
//    if (phred == "255")
//        return "NA";
//    else
//        return phred;
//}

std::vector<int> find_all_occurrence_in_string(std::string s) {
    std::vector<int> all_occurrence;
    std::string CIGAR_CHARACTER = "MSIDHNPX=";
    for (int i = 0; i < s.size(); i++) {
        size_t found = CIGAR_CHARACTER.find(s[i]);
        if (found != std::string::npos)
            all_occurrence.push_back(i);
    }
    return all_occurrence;
}

std::string find_location_or_cigar_on_read(int mapping_position, int query_position, std::string CIGAR_string, int max_length,
        bool isLocation, std::string BOTH_ADD_CIGAR, std::string REFERENCE_ADD_CIGAR, std::string READ_ADD_CIGAR) {

    int read_current_position = 0;
    int reference_current_position = mapping_position;
    std::vector<int> cigar_split = find_all_occurrence_in_string(CIGAR_string);
    cigar_split.insert(cigar_split.begin(), -1);
    std::string cigar_case;

    for (std::size_t i = 0; i < cigar_split.size()-1; i++) {
        int num_bases = std::stoi(CIGAR_string.substr(cigar_split[i]+1, cigar_split[i+1]-cigar_split[i]-1));
        cigar_case = CIGAR_string[cigar_split[i+1]];
        if (BOTH_ADD_CIGAR.find(cigar_case) != std::string::npos) {
            read_current_position += num_bases;
            reference_current_position += num_bases;
            if ((read_current_position > max_length-1) && (query_position > reference_current_position))
                return "NA";
            if (reference_current_position >= query_position) {
                if (read_current_position < reference_current_position - query_position)
                    return "NA";
                else {
                    int position_on_read = read_current_position - (reference_current_position - query_position);
                    if (position_on_read > max_length-1)
                        return "NA";
                    else {
                        if (isLocation)
                            return std::to_string(position_on_read);
                        else
                            return cigar_case;
                    }
                }
            }
        }
        else if (REFERENCE_ADD_CIGAR.find(cigar_case) != std::string::npos) {
            reference_current_position += num_bases;
            if (reference_current_position >= query_position) {
                if (isLocation)
                    return "D";
                else
                    return cigar_case;
            }
        }
        else if (READ_ADD_CIGAR.find(cigar_case) != std::string::npos) {
            read_current_position += num_bases;
            if (read_current_position > max_length - 1)
                return "NA";
        }
    }
    return "NA";
}

std::vector<std::string> extract_cigar_base(int mapping_position, int query_position, std::string CIGAR_string, int max_length, int window) {

    std::vector<std::string> base;
    for (int d = -window; d < window + 1; d++)
        base.push_back(find_location_or_cigar_on_read(mapping_position, query_position + d, CIGAR_string, max_length, false));
    return base;
}

std::string extract_base_quality_given_distance(int variant_position, int mapping_position, int distance, std::string quality_string, std::string CIGAR_string) {
    std::string location_on_read = find_location_or_cigar_on_read(mapping_position, variant_position + distance, CIGAR_string, quality_string.length(), true);
    if (location_on_read == "NA")
        return "NA";
    else if (location_on_read == "D")
        return "D";
    else {
        unsigned char quality_base = quality_string[std::stoi(location_on_read)];
        int ascii = quality_base;
        return std::to_string(ascii-33);
    }
}

std::vector<std::string> extract_base_quality(int variant_position, int mapping_position, std::string quality_string, std::string CIGAR_string, int window) {

    std::vector<std::string> quality;
    for (int d = -window; d < window + 1; d++)
        quality.push_back(extract_base_quality_given_distance(variant_position, mapping_position, d, quality_string, CIGAR_string));
    return quality;
}

std::string extract_reference_base_read_base_given_distance(int variant_position, int mapping_position, int distance,
        std::string base_string, std::string reference_string, std::string CIGAR_string, int window) {

    char ref = reference_string[distance + window];
    std::string location_on_read = find_location_or_cigar_on_read(mapping_position, variant_position + distance, CIGAR_string, base_string.length(), true);
    std::string result;
    if (location_on_read == "NA")
        return "NA";
    else if (location_on_read == "D") {
        result.push_back(ref);
        result.push_back('|');
        result.push_back('D');
        return result;
    } else {
        result.push_back(ref);
        result.push_back('|');
        result.push_back(base_string[std::stoi(location_on_read)]);
        return result;
    }
}

std::vector<std::string> extract_reference_base_read_base(int variant_position, int mapping_position, std::string base_string, std::string reference_string, std::string CIGAR_string, int window) {

    std::vector<std::string> bases;
    for (int d = -window; d < window + 1; d++)
        bases.push_back(extract_reference_base_read_base_given_distance(variant_position, mapping_position, d, base_string, reference_string, CIGAR_string, window));
    return bases;
}

std::vector<int> find_indel_position(std::vector<int> &cigar_split, std::string cigar_string, int position) {

    std::vector<int> insertion_start;
    int current = position;
    std::vector<int> string_index = cigar_split;
    string_index.insert(string_index.begin(),-1);
    std::string refer_string = "MX=N";

    for (int i = 0; i < cigar_split.size(); i++) {
        int id = cigar_split[i];
        int length = std::stoi(cigar_string.substr(string_index[i]+1, string_index[i+1]-string_index[i]-1));
        std::string base;
        base.push_back(cigar_string[id]);

        if (base == "I") {
            insertion_start.push_back(current - 1);
            insertion_start.push_back(current);
        }
        else if (base == "D") {
            int start = current - 1;
            current += length;
            int end = current + 1;
            for (int j = start; j < end; j++)
                insertion_start.push_back(j);
        }
        else if (refer_string.find(base) != std::string::npos)
            current += length;
    }
    return insertion_start;
}

std::string extract_nearby_indel(int variant_position, int mapping_position, std::string CIGAR_string) {

    std::vector<int> cigar_split = find_all_occurrence_in_string(CIGAR_string);
    std::vector<int> all_location_involved_in_indel = find_indel_position(cigar_split, CIGAR_string, mapping_position);

    if (all_location_involved_in_indel.empty())
        return "-1";
    else {
        std::vector<int> distance_to_variant_position;
        for (int &l : all_location_involved_in_indel)
            distance_to_variant_position.push_back(std::abs(l - variant_position));
        return std::to_string(*std::min_element(distance_to_variant_position.begin(), distance_to_variant_position.end()));

    }
}

std::string extract_homopolymer_on_read(int variant_position, int mapping_position, std::string base_string, std::string CIGAR_string,
        int HOMOPOLYMER_SIZE, const std::vector<std::string> &HOMOPOLYMER) {

    std::string location_on_read = find_location_or_cigar_on_read(mapping_position, variant_position, CIGAR_string, base_string.length(), true);
    if ((location_on_read == "D") || (location_on_read == "NA"))
        return "NA";
    std::string neighborhood;
    for (int i = - HOMOPOLYMER_SIZE + 1; i < HOMOPOLYMER_SIZE; i++) {
        location_on_read = find_location_or_cigar_on_read(mapping_position, variant_position + i, CIGAR_string, base_string.length(), true);
        if ((location_on_read == "NA") || (location_on_read == "D"))
            neighborhood.push_back('0');
        else
            neighborhood.push_back(base_string[std::stoi(location_on_read)]);
    }
    for (int i = 0; i < HOMOPOLYMER_SIZE; i++) {
        std::string sub_string = neighborhood.substr(i, HOMOPOLYMER_SIZE);
        if (std::find(HOMOPOLYMER.begin(), HOMOPOLYMER.end(), sub_string) != HOMOPOLYMER.end())
            return "1";
    }
    return "0";
}

std::string extract_homopolymer_on_reference(std::string ref_string, int HOMOPOLYMER_SIZE, const std::vector<std::string> &HOMOPOLYMER) {

    for (int i = 0; i < ref_string.length()-HOMOPOLYMER_SIZE+1; i++) {
        std::string sub_string = ref_string.substr(i, HOMOPOLYMER_SIZE);
        if (std::find(HOMOPOLYMER.begin(), HOMOPOLYMER.end(), sub_string) != HOMOPOLYMER.end()) {
            return "1";
        }

    }
    return "0";
}

//std::string extract_CIGAR_string(std::string read) {
//
//    return split(read, "\t")[5];
//}

std::string replaceChar(std::string str, char ch1, char ch2) {

    for (int i = 0; i < str.length(); ++i) {
        if (str[i] == ch1) {
            str[i] = ch2;
            break; // Only one x?
        }

    }
    return str;
}

void read_preparebed(const std::vector<std::vector<std::string>> &preparebed, std::map<std::string, std::vector<int>> &variant_position_dict,
        std::map<std::string, std::vector<std::string>> &variant_base_dict, std::map<std::string, std::vector<std::string>> &context_dict) {

    std::string chr, pos, ctx, ref, var;
    for (auto const &sp : preparebed) {
        chr = sp[0];
        pos = sp[2];
        ctx = sp[3];
        ref = sp[4];
        var = sp[5];

        variant_position_dict[chr].push_back(std::stoi(pos));
        variant_base_dict[chr].push_back(var);

//        char ref_char[1];
//        std::strcpy(ref_char, ref.c_str());
//        std::replace(ctx.begin(), ctx.end(), 'x', ref_char[0]);
        char ref_char = ref[0];
        std::string ctx_new = replaceChar(ctx, 'x', ref_char);
        context_dict[chr].push_back(ctx_new);

    }
}

void extract_features_from_reads_filter_cluster(std::string sam_file, std::string out_paired_reads_qsort_features,
                                                std::map<std::string, std::vector<int>> &variant_position_dict,
                                                std::map<std::string, std::vector<std::string>> &variant_base_dict,
                                                std::map<std::string, std::vector<std::string>> &context_dict,
                                                std::map<std::string, int> &mapping_to_var,
                                                std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<int>>>> &matrix,
                                                std::map<int, std::vector<int>> &read_pos_loc,
                                                int HOMOPOLYMER_SIZE, const std::vector<std::string> &HOMOPOLYMER, int window) {

    std::ofstream outfile(out_paired_reads_qsort_features);
    outfile.close();
    outfile.open(out_paired_reads_qsort_features, std::ios_base::app);

    std::ifstream inFile;
    std::string tmp0, tmp1;
    std::vector<std::string> sp0, sp1;
    std::string chr1, chr0;
    std::string pos1, pos0;
    std::string basestring1, basestring0;
    std::string qualstring1, qualstring0;
    std::string CIGAR1, CIGAR0;
    std::string flag1, flag0;
    std::string binflag1, binflag0;
    std::string strand1, strand0;
    std::string MQ1, MQ0;
    std::vector<std::string> INFO_1, INFO_0;
    std::vector<std::string> SUP_INFO_1, SUP_INFO_0;
    int size, samechr;
    std::vector<std::string> RELATION;
    std::vector<int> related_variant, related_read, related_position;
    int matrix_id, row_id, col_id;

    inFile.open(sam_file);
    std::getline(inFile, tmp0);
    sp0 = split(tmp0, "\t");
    while (!tmp0.empty()) {
        tmp1 = tmp0;
        sp1 = sp0;
        std::getline(inFile, tmp0);
        if (tmp0.empty())
            break;
        sp0 = split(tmp0, "\t");
        if (sp0[0] != sp1[0])
            continue;
        chr1 = sp1[2];
        chr0 = sp0[2];
        pos1 = sp1[3];
        pos0 = sp0[3];
        basestring1 = sp1[9];
        basestring0 = sp0[9];
        qualstring1 = sp1[10];
        qualstring0 = sp0[10];
        CIGAR1 = sp1[5];
        CIGAR0 = sp0[5];
        if ((CIGAR1 == "*") || (CIGAR0 == "*"))
            continue;
        flag1 = sp1[1];
        flag0 = sp0[1];
        binflag1 = std::bitset<8>(std::stoi(flag1)).to_string();
        binflag0 = std::bitset<8>(std::stoi(flag0)).to_string();
        strand1 = binflag1[3];
        strand0 = binflag0[3];
        MQ1 = sp1[4];
        MQ0 = sp0[4];
        INFO_1.assign({MQ1, CIGAR1});
        INFO_0.assign({MQ0, CIGAR0});
        SUP_INFO_1.assign({chr1, pos1});
        SUP_INFO_0.assign({chr0, pos0});


        if (chr1 == chr0) {
            related_variant.clear();
            related_read.clear();
            related_position.clear();

            size = std::abs(std::stoi(pos0) - std::stoi(pos1)) + basestring1.length();
            samechr = 1;
            RELATION.assign({std::to_string(size), std::to_string(samechr)});
            for (int i_var = 0; i_var < variant_position_dict[chr1].size(); i_var++) {
                int var_pos = variant_position_dict[chr1][i_var];
                std::string var_loc_on_read1(find_location_or_cigar_on_read(std::stoi(pos1), var_pos, CIGAR1, basestring1.length(), true));
                if (var_loc_on_read1 != "NA") {
                    if (var_loc_on_read1 == "D")
                        continue;
                    std::string var_base1;
                    var_base1 += basestring1[std::stoi(var_loc_on_read1)];
                    if (var_base1 == variant_base_dict[chr1][i_var]) {
                        related_variant.push_back(i_var);
                        related_read.push_back(1);
                        related_position.push_back(std::stoi(var_loc_on_read1));
                        continue;
                    }
                }
                std::string var_loc_on_read0(find_location_or_cigar_on_read(std::stoi(pos0), var_pos, CIGAR0, basestring0.length(), true));

                if (var_loc_on_read0 != "NA") {
                    if (var_loc_on_read0 == "D")
                        continue;
                    std::string var_base0;
                    var_base0 += basestring0[std::stoi(var_loc_on_read0)];
                    if (var_base0 == variant_base_dict[chr0][i_var]) {
                        related_variant.push_back(i_var);
                        related_read.push_back(0);
                        related_position.push_back(std::stoi(var_loc_on_read0));
                        read_pos_loc[var_pos].push_back(std::stoi(var_loc_on_read0));
                        continue;
                    }
                }
            }
            for (int &rel_var1 : related_variant) {
                std::string rel_var1_string(chr1 + "-" + std::to_string(variant_position_dict[chr1][rel_var1]));
                matrix_id = mapping_to_var[rel_var1_string];
                std::vector<std::string> matrix_i;

                if (matrix_id != -1) {
                    matrix_i = matrix[matrix_id].first;
                    auto it1 = std::find(matrix_i.begin(), matrix_i.end(), rel_var1_string);
                    if (it1 != matrix_i.end())
                        row_id = it1 - matrix_i.begin();
                    else
                        continue;
                    matrix[matrix_id].second[row_id][row_id] += 1;
                } else
                    continue;

                for (int &rel_var2 : related_variant) {
                    std::string rel_var2_string(chr1 + "-" + std::to_string(variant_position_dict[chr1][rel_var2]));
                    if ((mapping_to_var[rel_var2_string] != -1) && (rel_var2_string != rel_var1_string)) {
                        auto it2 = std::find(matrix_i.begin(), matrix_i.end(), rel_var2_string);
                        if (it2 != matrix_i.end())
                            col_id = it2 - matrix_i.begin();
                        else
                            continue;
                        matrix[matrix_id].second[row_id][col_id] += 1;
                    } else
                        continue;
                }
            }
            for (int i_rel_var = 0; i_rel_var < related_variant.size(); i_rel_var++) {
                int i_var = related_variant[i_rel_var];
                int variant_position = variant_position_dict[chr1][i_var];
                std::string reference_string(context_dict[chr1][i_var]);
                std::vector<std::string> INFOA_1, INFOA_0;

                std::vector<std::string> VARIANT = {chr1, std::to_string(variant_position), reference_string};
                std::vector<std::string> EXTRA = {extract_homopolymer_on_reference(reference_string, HOMOPOLYMER_SIZE, HOMOPOLYMER)};
                std::vector<std::string> quality1 = extract_base_quality(variant_position, std::stoi(pos1), qualstring1, CIGAR1, window);
                std::vector<std::string> quality0 = extract_base_quality(variant_position, std::stoi(pos0), qualstring0, CIGAR0, window);
                std::vector<std::string> base1 = extract_reference_base_read_base(variant_position, std::stoi(pos1), basestring1, reference_string, CIGAR1, window);
                std::vector<std::string> base0 = extract_reference_base_read_base(variant_position, std::stoi(pos0), basestring0, reference_string, CIGAR0, window);
                std::vector<std::string> cigar1 = extract_cigar_base(std::stoi(pos1), variant_position, CIGAR1, basestring1.length(), window);
                std::vector<std::string> cigar0 = extract_cigar_base(std::stoi(pos0), variant_position, CIGAR0, basestring0.length(), window);
                std::string indel1 = extract_nearby_indel(variant_position, std::stoi(pos1), CIGAR1);
                std::string indel0 = extract_nearby_indel(variant_position, std::stoi(pos0), CIGAR0);
                std::string homopolymer1 = extract_homopolymer_on_read(variant_position, std::stoi(pos1), basestring1, CIGAR1, HOMOPOLYMER_SIZE, HOMOPOLYMER);
                std::string homopolymer0 = extract_homopolymer_on_read(variant_position, std::stoi(pos0), basestring0, CIGAR0, HOMOPOLYMER_SIZE, HOMOPOLYMER);
                std::vector<std::string> indel_homo1 = {indel1, homopolymer1};
                std::vector<std::string> indel_homo0 = {indel0, homopolymer0};

                INFOA_1.reserve(INFO_1.size() + quality1.size() + base1.size() + cigar1.size() + indel_homo1.size());
                INFOA_1.insert(std::end(INFOA_1), std::begin(INFO_1), std::end(INFO_1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(quality1), std::end(quality1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(base1), std::end(base1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(cigar1), std::end(cigar1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(indel_homo1), std::end(indel_homo1));

                INFOA_0.reserve(INFO_0.size() + quality0.size() + base0.size() + cigar0.size() + indel_homo0.size());
                INFOA_0.insert(std::end(INFOA_0), std::begin(INFO_0), std::end(INFO_0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(quality0), std::end(quality0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(base0), std::end(base0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(cigar0), std::end(cigar0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(indel_homo0), std::end(indel_homo0));

                int related_pos = related_position[i_rel_var];

                if (related_read[i_rel_var] == 0) {
                    std::vector<std::string> RELATIONA, A = {strand0, strand1, flag0, flag1}, out_list, rest1;
                    RELATIONA.reserve(RELATION.size() + A.size());
                    RELATIONA.insert(std::end(RELATIONA), std::begin(RELATION), std::end(RELATION));
                    RELATIONA.insert(std::end(RELATIONA), std::begin(A), std::end(A));
                    size_t find_h = CIGAR0.find('H');
                    size_t find_s = CIGAR0.find('S');
                    std::string find_s_str(std::to_string(find_s));
                    if (find_s == std::string::npos)
                        find_s_str = "-1";

                    rest1 = {std::to_string(related_pos), std::to_string(basestring0.length()-related_pos), find_s_str};
                    out_list.reserve(INFOA_0.size() + RELATIONA.size() + INFOA_1.size() + EXTRA.size() + SUP_INFO_0.size() + SUP_INFO_1.size() + VARIANT.size() + rest1.size());
                    out_list.insert(std::end(out_list), std::begin(INFOA_0), std::end(INFOA_0));
                    out_list.insert(std::end(out_list), std::begin(RELATIONA), std::end(RELATIONA));
                    out_list.insert(std::end(out_list), std::begin(INFOA_1), std::end(INFOA_1));
                    out_list.insert(std::end(out_list), std::begin(EXTRA), std::end(EXTRA));
                    out_list.insert(std::end(out_list), std::begin(SUP_INFO_0), std::end(SUP_INFO_0));
                    out_list.insert(std::end(out_list), std::begin(SUP_INFO_1), std::end(SUP_INFO_1));
                    out_list.insert(std::end(out_list), std::begin(VARIANT), std::end(VARIANT));
                    out_list.insert(std::end(out_list), std::begin(rest1), std::end(rest1));

                    if ((find_s_str != "-1") && (find_s_str != "0"))
                        continue;
                    if ((find_h != std::string::npos) && (find_h != 0))
                        continue;

                    for (std::string &e : out_list)
                        outfile << e << "\t";
                    outfile << "\n";
                }
                else {
                    std::vector<std::string> RELATIONA, B = {strand1, strand0, flag1, flag0}, out_list, rest2;
                    RELATIONA.reserve(RELATION.size() + B.size());
                    RELATIONA.insert(std::end(RELATIONA), std::begin(RELATION), std::end(RELATION));
                    RELATIONA.insert(std::end(RELATIONA), std::begin(B), std::end(B));
                    size_t find_h = CIGAR1.find('H');
                    size_t find_s = CIGAR1.find('S');
                    std::string find_s_str(std::to_string(find_s));
                    if (find_s == std::string::npos)
                        find_s_str = "-1";

                    rest2 = {std::to_string(related_pos), std::to_string(basestring1.length()-related_pos), find_s_str};
                    out_list.reserve(INFOA_1.size() + RELATIONA.size() + INFOA_0.size() + EXTRA.size() + SUP_INFO_1.size() + SUP_INFO_0.size() + VARIANT.size() + rest2.size());
                    out_list.insert(std::end(out_list), std::begin(INFOA_1), std::end(INFOA_1));
                    out_list.insert(std::end(out_list), std::begin(RELATIONA), std::end(RELATIONA));
                    out_list.insert(std::end(out_list), std::begin(INFOA_0), std::end(INFOA_0));
                    out_list.insert(std::end(out_list), std::begin(EXTRA), std::end(EXTRA));
                    out_list.insert(std::end(out_list), std::begin(SUP_INFO_1), std::end(SUP_INFO_1));
                    out_list.insert(std::end(out_list), std::begin(SUP_INFO_0), std::end(SUP_INFO_0));
                    out_list.insert(std::end(out_list), std::begin(VARIANT), std::end(VARIANT));
                    out_list.insert(std::end(out_list), std::begin(rest2), std::end(rest2));

                    if ((find_s_str != "-1") && (find_s_str != "0"))
                        continue;
                    if ((find_h != std::string::npos) && (find_h != 0))
                        continue;

                    for (std::string &e : out_list)
                        outfile << e << "\t";
                    outfile << "\n";

                }
            }
        }
        else {
            size = -1;
            samechr = 0;
            RELATION.assign({std::to_string(size), std::to_string(samechr)});
            std::vector<int> related_variant0, related_variant1, related_position0, related_position1;
            related_variant0.clear();
            related_variant1.clear();
            related_position0.clear();
            related_position1.clear();

            for (int i_var = 0; i_var < variant_position_dict[chr1].size(); i_var++) {
                int var_pos = variant_position_dict[chr1][i_var];
                std::string var_loc_on_read1(find_location_or_cigar_on_read(std::stoi(pos1), var_pos, CIGAR1, basestring1.length(), true));
                if (var_loc_on_read1 != "NA") {
                    if (var_loc_on_read1 == "D")
                        continue;
                    std::string var_base1;
                    var_base1 += basestring1[std::stoi(var_loc_on_read1)];
                    if (var_base1 == variant_base_dict[chr1][i_var]) {
                        related_variant1.push_back(i_var);
                        related_position1.push_back(std::stoi(var_loc_on_read1));
                        read_pos_loc[var_pos].push_back(std::stoi(var_loc_on_read1));
                        continue;
                    }
                }
            }
            for (int &rel_var1 : related_variant1) {
                std::string rel_var1_string(chr1 + "-" + std::to_string(variant_position_dict[chr1][rel_var1]));
                matrix_id = mapping_to_var[rel_var1_string];
                std::vector<std::string> matrix_i;

                if (matrix_id != -1) {
                    matrix_i = matrix[matrix_id].first;
                    auto it1 = std::find(matrix_i.begin(), matrix_i.end(), rel_var1_string);
                    if (it1 != matrix_i.end())
                        row_id = it1 - matrix_i.begin();
                    else
                        continue;
                    matrix[matrix_id].second[row_id][row_id] += 1;
                } else
                    continue;
                for (int &rel_var2 : related_variant1) {
                    std::string rel_var2_string(chr1 + "-" + std::to_string(variant_position_dict[chr1][rel_var2]));
                    if ((mapping_to_var[rel_var2_string] != -1) && (rel_var2_string != rel_var1_string)) {
                        auto it2 = std::find(matrix_i.begin(), matrix_i.end(), rel_var2_string);
                        if (it2 != matrix_i.end())
                            col_id = it2 - matrix_i.begin();
                        else
                            continue;
                        matrix[matrix_id].second[row_id][col_id] += 1;
                    } else
                        continue;
                }
            }
            for (int i_var = 0; i_var < variant_position_dict[chr0].size(); i_var++) {
                int var_pos = variant_position_dict[chr0][i_var];
                std::string var_loc_on_read0(find_location_or_cigar_on_read(std::stoi(pos0), var_pos, CIGAR0, basestring0.length(), true));
                if (var_loc_on_read0 != "NA") {
                    if (var_loc_on_read0 == "D")
                        continue;
                    std::string var_base0;
                    var_base0 += basestring0[std::stoi(var_loc_on_read0)];
                    if (var_base0 == variant_base_dict[chr0][i_var]) {
                        related_variant0.push_back(i_var);
                        related_position0.push_back(std::stoi(var_loc_on_read0));
                        continue;
                    }
                }
            }
            for (int &rel_var1 : related_variant0) {
                std::string rel_var1_string(chr0 + "-" + std::to_string(variant_position_dict[chr0][rel_var1]));
                matrix_id = mapping_to_var[rel_var1_string];
                std::vector<std::string> matrix_i;

                if (matrix_id != -1) {
                    matrix_i = matrix[matrix_id].first;
                    auto it1 = std::find(matrix_i.begin(), matrix_i.end(), rel_var1_string);
                    if (it1 != matrix_i.end())
                        row_id = it1 - matrix_i.begin();
                    else
                        continue;
                    matrix[matrix_id].second[row_id][row_id] += 1;
                } else
                    continue;
                for (int &rel_var2 : related_variant0) {
                    std::string rel_var2_string(chr0 + "-" + std::to_string(variant_position_dict[chr0][rel_var2]));
                    if ((mapping_to_var[rel_var2_string] != -1) && (rel_var2_string != rel_var1_string)) {
                        auto it2 = std::find(matrix_i.begin(), matrix_i.end(), rel_var2_string);
                        if (it2 != matrix_i.end())
                            col_id = it2 - matrix_i.begin();
                        else
                            continue;
                        matrix[matrix_id].second[row_id][col_id] += 1;
                    } else
                        continue;
                }
            }
            for (int i_rel_var = 0; i_rel_var < related_variant1.size(); i_rel_var++) {
                int i_var = related_variant1[i_rel_var];
                int variant_position = variant_position_dict[chr1][i_var];
                std::string reference_string(context_dict[chr1][i_var]);
                std::vector<std::string> INFOA_1, INFOA_0;

                std::vector<std::string> VARIANT = {chr1, std::to_string(variant_position), reference_string};
                std::vector<std::string> EXTRA = {extract_homopolymer_on_reference(reference_string, HOMOPOLYMER_SIZE, HOMOPOLYMER)};
                std::vector<std::string> quality1 = extract_base_quality(variant_position, std::stoi(pos1), qualstring1, CIGAR1, window);
                std::vector<std::string> quality0 = extract_base_quality(variant_position, std::stoi(pos0), qualstring0, CIGAR0, window);
                std::vector<std::string> base1 = extract_reference_base_read_base(variant_position, std::stoi(pos1), basestring1, reference_string, CIGAR1, window);
                std::vector<std::string> base0 = extract_reference_base_read_base(variant_position, std::stoi(pos0), basestring0, reference_string, CIGAR0, window);
                std::vector<std::string> cigar1 = extract_cigar_base(std::stoi(pos1), variant_position, CIGAR1, basestring1.length(), window);
                std::vector<std::string> cigar0 = extract_cigar_base(std::stoi(pos0), variant_position, CIGAR0, basestring0.length(), window);
                std::string indel1 = extract_nearby_indel(variant_position, std::stoi(pos1), CIGAR1);
                std::string indel0 = extract_nearby_indel(variant_position, std::stoi(pos0), CIGAR0);
                std::string homopolymer1 = extract_homopolymer_on_read(variant_position, std::stoi(pos1), basestring1, CIGAR1, HOMOPOLYMER_SIZE, HOMOPOLYMER);
                std::string homopolymer0 = extract_homopolymer_on_read(variant_position, std::stoi(pos0), basestring0, CIGAR0, HOMOPOLYMER_SIZE, HOMOPOLYMER);
                std::vector<std::string> indel_homo1 = {indel1, homopolymer1};
                std::vector<std::string> indel_homo0 = {indel0, homopolymer0};

                INFOA_1.reserve(INFO_1.size() + quality1.size() + base1.size() + cigar1.size() + indel_homo1.size());
                INFOA_1.insert(std::end(INFOA_1), std::begin(INFO_1), std::end(INFO_1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(quality1), std::end(quality1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(base1), std::end(base1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(cigar1), std::end(cigar1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(indel_homo1), std::end(indel_homo1));

                INFOA_0.reserve(INFO_0.size() + quality0.size() + base0.size() + cigar0.size() + indel_homo0.size());
                INFOA_0.insert(std::end(INFOA_0), std::begin(INFO_0), std::end(INFO_0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(quality0), std::end(quality0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(base0), std::end(base0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(cigar0), std::end(cigar0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(indel_homo0), std::end(indel_homo0));

                std::vector<std::string> RELATIONA, out_list;
                std::vector<std::string> B = {strand1, strand0, flag1, flag0};
                RELATIONA.reserve(RELATION.size() + B.size());
                RELATIONA.insert(std::end(RELATIONA), std::begin(RELATION), std::end(RELATION));
                RELATIONA.insert(std::end(RELATIONA), std::begin(B), std::end(B));
                size_t find_h = CIGAR1.find('H');
                size_t find_s = CIGAR1.find('S');
                std::string find_s_str = std::to_string(find_s);
                if (find_s == std::string::npos)
                    find_s_str = "-1";

                int related_pos = related_position1[i_rel_var];
                std::vector<std::string> rest1 = {std::to_string(related_pos), std::to_string(basestring1.length()-related_pos), find_s_str};
                out_list.reserve(INFOA_1.size() + RELATIONA.size() + INFOA_0.size() + EXTRA.size() + SUP_INFO_1.size() + SUP_INFO_0.size() + VARIANT.size() + rest1.size());
                out_list.insert(std::end(out_list), std::begin(INFOA_1), std::end(INFOA_1));
                out_list.insert(std::end(out_list), std::begin(RELATIONA), std::end(RELATIONA));
                out_list.insert(std::end(out_list), std::begin(INFOA_0), std::end(INFOA_0));
                out_list.insert(std::end(out_list), std::begin(EXTRA), std::end(EXTRA));
                out_list.insert(std::end(out_list), std::begin(SUP_INFO_1), std::end(SUP_INFO_1));
                out_list.insert(std::end(out_list), std::begin(SUP_INFO_0), std::end(SUP_INFO_0));
                out_list.insert(std::end(out_list), std::begin(VARIANT), std::end(VARIANT));
                out_list.insert(std::end(out_list), std::begin(rest1), std::end(rest1));

                if ((find_s_str != "-1") && (find_s_str != "0"))
                    continue;
                if ((find_h != std::string::npos) && (find_h != 0))
                    continue;

                for (std::string &e : out_list)
                    outfile << e << "\t";
                outfile << "\n";
            }
            for (int i_rel_var = 0; i_rel_var < related_variant0.size(); i_rel_var++) {
                int i_var = related_variant0[i_rel_var];
                int variant_position = variant_position_dict[chr0][i_var];
                std::string reference_string(context_dict[chr0][i_var]);
                std::vector<std::string> INFOA_1, INFOA_0;

                std::vector<std::string> VARIANT = {chr0, std::to_string(variant_position), reference_string};
                std::vector<std::string> EXTRA = {extract_homopolymer_on_reference(reference_string, HOMOPOLYMER_SIZE, HOMOPOLYMER)};
                std::vector<std::string> quality1 = extract_base_quality(variant_position, std::stoi(pos1), qualstring1, CIGAR1, window);
                std::vector<std::string> quality0 = extract_base_quality(variant_position, std::stoi(pos0), qualstring0, CIGAR0, window);
                std::vector<std::string> base1 = extract_reference_base_read_base(variant_position, std::stoi(pos1), basestring1, reference_string, CIGAR1, window);
                std::vector<std::string> base0 = extract_reference_base_read_base(variant_position, std::stoi(pos0), basestring0, reference_string, CIGAR0, window);
                std::vector<std::string> cigar1 = extract_cigar_base(std::stoi(pos1), variant_position, CIGAR1, basestring1.length(), window);
                std::vector<std::string> cigar0 = extract_cigar_base(std::stoi(pos0), variant_position, CIGAR0, basestring0.length(), window);
                std::string indel1 = extract_nearby_indel(variant_position, std::stoi(pos1), CIGAR1);
                std::string indel0 = extract_nearby_indel(variant_position, std::stoi(pos0), CIGAR0);
                std::string homopolymer1 = extract_homopolymer_on_read(variant_position, std::stoi(pos1), basestring1, CIGAR1, HOMOPOLYMER_SIZE, HOMOPOLYMER);
                std::string homopolymer0 = extract_homopolymer_on_read(variant_position, std::stoi(pos0), basestring0, CIGAR0, HOMOPOLYMER_SIZE, HOMOPOLYMER);
                std::vector<std::string> indel_homo1 = {indel1, homopolymer1};
                std::vector<std::string> indel_homo0 = {indel0, homopolymer0};

                INFOA_1.reserve(INFO_1.size() + quality1.size() + base1.size() + cigar1.size() + indel_homo1.size());
                INFOA_1.insert(std::end(INFOA_1), std::begin(INFO_1), std::end(INFO_1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(quality1), std::end(quality1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(base1), std::end(base1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(cigar1), std::end(cigar1));
                INFOA_1.insert(std::end(INFOA_1), std::begin(indel_homo1), std::end(indel_homo1));

                INFOA_0.reserve(INFO_0.size() + quality0.size() + base0.size() + cigar0.size() + indel_homo0.size());
                INFOA_0.insert(std::end(INFOA_0), std::begin(INFO_0), std::end(INFO_0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(quality0), std::end(quality0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(base0), std::end(base0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(cigar0), std::end(cigar0));
                INFOA_0.insert(std::end(INFOA_0), std::begin(indel_homo0), std::end(indel_homo0));

                std::vector<std::string> RELATIONA, out_list;
                std::vector<std::string> A = {strand0, strand1, flag0, flag1};
                RELATIONA.reserve(RELATION.size() + A.size());
                RELATIONA.insert(std::end(RELATIONA), std::begin(RELATION), std::end(RELATION));
                RELATIONA.insert(std::end(RELATIONA), std::begin(A), std::end(A));
                size_t find_h = CIGAR0.find('H');
                size_t find_s = CIGAR0.find('S');
                std::string find_s_str = std::to_string(find_s);
                if (find_s == std::string::npos)
                    find_s_str = "-1";

                int related_pos = related_position0[i_rel_var];
                std::vector<std::string> rest2 = {std::to_string(related_pos), std::to_string(basestring0.length()-related_pos), find_s_str};
                out_list.reserve(INFOA_0.size() + RELATIONA.size() + INFOA_1.size() + EXTRA.size() + SUP_INFO_0.size() + SUP_INFO_1.size() + VARIANT.size() + rest2.size());
                out_list.insert(std::end(out_list), std::begin(INFOA_0), std::end(INFOA_0));
                out_list.insert(std::end(out_list), std::begin(RELATIONA), std::end(RELATIONA));
                out_list.insert(std::end(out_list), std::begin(INFOA_1), std::end(INFOA_1));
                out_list.insert(std::end(out_list), std::begin(EXTRA), std::end(EXTRA));
                out_list.insert(std::end(out_list), std::begin(SUP_INFO_0), std::end(SUP_INFO_0));
                out_list.insert(std::end(out_list), std::begin(SUP_INFO_1), std::end(SUP_INFO_1));
                out_list.insert(std::end(out_list), std::begin(VARIANT), std::end(VARIANT));
                out_list.insert(std::end(out_list), std::begin(rest2), std::end(rest2));

                if ((find_s_str != "-1") && (find_s_str != "0"))
                    continue;
                if ((find_h != std::string::npos) && (find_h != 0))
                    continue;

                for (std::string &e : out_list)
                    outfile << e << "\t";
                outfile << "\n";
            }

        }
    }
    outfile.close();
    inFile.close();
}

void write_to_filter_cluster_bed(std::string filter_cluster_bed, std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<int>>>> &matrix,
                                 std::map<std::string, std::vector<int>> &variant_position_dict, std::map<int, std::vector<int>> &read_pos_loc,
                                 double CLUSTERED_FRAC, int NONCLUSTERED_COUNT) {

    std::ofstream out_var_txt(filter_cluster_bed+".var.txt");
    out_var_txt.close();
    out_var_txt.open(filter_cluster_bed+".var.txt", std::ios_base::app);
    std::vector<std::string> clustered_variants;

    for (std::pair<std::vector<std::string>, std::vector<std::vector<int>>> &item : matrix) {
        for (std::string &e : item.first)
            out_var_txt << e << "\t";
        out_var_txt << "\n";

        for (int i_var_main = 0; i_var_main < item.first.size(); i_var_main++) {
            int cnt_var_main = item.second[i_var_main][i_var_main];
            if (cnt_var_main == 0)
                continue;
            int flag = 0;
            int pos = std::stoi(split(item.first[i_var_main], "-")[1]);

            for (int i_var_else = 0; i_var_else < item.first.size(); i_var_else++) {
                if (i_var_main == i_var_else)
                    continue;
                int pos_else = std::stoi(split(item.first[i_var_else], "-")[1]);
                if (std::abs(pos_else - pos) < 3)
                    continue;
                int cnt_cooccur = item.second[i_var_main][i_var_else];
                if (((cnt_cooccur*1.0)/(cnt_var_main*1.0) > CLUSTERED_FRAC) && ((cnt_var_main-cnt_cooccur) < NONCLUSTERED_COUNT)) {
                    flag = 1;
                }
            }
            if (flag == 1)
                clustered_variants.push_back(item.first[i_var_main]);
        }
    }
    out_var_txt.close();

    std::ofstream out_bed(filter_cluster_bed);
    out_bed.close();
    out_bed.open(filter_cluster_bed, std::ios_base::app);

    for (std::string &var : clustered_variants) {
        std::vector<std::string> sp = split(var, "-");
        sp.push_back(sp[1]);
        sp.at(1) = std::to_string(std::stoi(sp[1])-1);

        for (int i = 0; i < 3; i++) {
            out_bed << sp[i];
            if (i != 2)
                out_bed << "\t";
        }
        out_bed << "\n";

        // for (std::string &e : sp)
        //     out_bed << e << "\t";
    }
    out_bed.close();
}

// [[Rcpp::export]]
int machinelearn_extract_features_from_reads_filter_cluster(std::string inpath, std::string outpath, std::string sample_id, int WINDOW) {

    //std::cout << "Begin py78" << std::endl;
    // std::string inpath = argv[1];
    // std::string outpath = argv[2];
    // std::string sample_id = argv[3];

    // int WINDOW = 3;
    std::vector<std::vector<std::string>> preparebed = read_bed_for_feature(inpath+"/"+sample_id+".fastabed", WINDOW);

    std::map<std::string, int> mapping_to_var;
    std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<int>>>> matrix;
    build_mapping_to_var(preparebed, mapping_to_var, matrix);

    int HOMOPOLYMER_SIZE = 5;
    std::vector<std::string> HOMOPOLYMER = make_homopolymer(HOMOPOLYMER_SIZE);

    std::map<std::string, std::vector<int>> variant_position_dict;
    std::map<std::string, std::vector<std::string>> variant_base_dict;
    std::map<std::string, std::vector<std::string>> context_dict;

    read_preparebed(preparebed, variant_position_dict, variant_base_dict, context_dict);

    std::map<int, std::vector<int>> read_pos_loc;
    extract_features_from_reads_filter_cluster(inpath+"/"+sample_id+".paired-reads.qsort.sam", outpath+"/"+sample_id+".paired-reads.qsort.features",
                                               variant_position_dict, variant_base_dict, context_dict, mapping_to_var, matrix, read_pos_loc,
                                               HOMOPOLYMER_SIZE, HOMOPOLYMER, WINDOW);

    double CLUSTERED_FRAC = 0.9;
    int NONCLUSTERED_COUNT = 2;

    write_to_filter_cluster_bed(outpath+"/"+sample_id+".filter_cluster.bed", matrix, variant_position_dict, read_pos_loc, CLUSTERED_FRAC, NONCLUSTERED_COUNT);

    return 0;
}
