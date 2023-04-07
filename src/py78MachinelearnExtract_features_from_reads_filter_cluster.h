//
// Created by Ran Hu on 8/1/21.
//

#ifndef INC_78_PY78MACHINELEARNEXTRACT_FEATURES_FROM_READS_FILTER_CLUSTER_H
#define INC_78_PY78MACHINELEARNEXTRACT_FEATURES_FROM_READS_FILTER_CLUSTER_H

#include <Rcpp.h>
using namespace Rcpp;

#include <map>

std::vector<std::string> make_homopolymer(int HOMOPOLYMER_SIZE);
std::vector<std::vector<std::string>> read_bed_for_feature(std::string fastabed, int window);
std::string trans_to_upper(std::string s);
void build_mapping_to_var(const std::vector<std::vector<std::string>> &preparebed, std::map<std::string, int> &mapping_to_var,
                          std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<int>>>> &matrix);

//std::string extract_mapping_quality(std::string read);
std::vector<int> find_all_occurrence_in_string(std::string s);

std::string find_location_or_cigar_on_read(int mapping_position, int query_position, std::string CIGAR_string, int max_length, bool isLocation,
                                           std::string BOTH_ADD_CIGAR="M=X", std::string REFERENCE_ADD_CIGAR="DN", std::string READ_ADD_CIGAR="IS");

std::vector<std::string> extract_cigar_base(int mapping_position, int query_position, std::string CIGAR_string, int max_length, int window);

std::string extract_base_quality_given_distance(int variant_position, int mapping_position, int distance, std::string quality_string, std::string CIGAR_string);
std::vector<std::string> extract_base_quality(int variant_position, int mapping_position, std::string quality_string, std::string CIGAR_string, int window);
std::string extract_reference_base_read_base_given_distance(int variant_position, int mapping_position, int distance,
                                                            std::string base_string, std::string reference_string, std::string CIGAR_string, int window);
std::vector<std::string> extract_reference_base_read_base(int variant_position, int mapping_position, std::string base_string, std::string reference_string, std::string CIGAR_string, int window);

std::vector<int> find_indel_position(std::vector<int> &cigar_split, std::string cigar_string, int position);
std::string extract_nearby_indel(int variant_position, int mapping_position, std::string CIGAR_string);
std::string extract_homopolymer_on_read(int variant_position, int mapping_position, std::string base_string, std::string CIGAR_string,
                                        int HOMOPOLYMER_SIZE, const std::vector<std::string> &HOMOPOLYMER);
std::string extract_homopolymer_on_reference(std::string ref_string, int HOMOPOLYMER_SIZE, const std::vector<std::string> &HOMOPOLYMER);
//std::string extract_CIGAR_string(std::string read);

std::string replaceChar(std::string str, char ch1, char ch2);

void read_preparebed(const std::vector<std::vector<std::string>> &preparebed, std::map<std::string, std::vector<int>> &variant_position_dict,
                     std::map<std::string, std::vector<std::string>> &variant_base_dict, std::map<std::string, std::vector<std::string>> &context_dict);

void extract_features_from_reads_filter_cluster(std::string sam_file, std::string out_paired_reads_qsort_features,
                                                std::map<std::string, std::vector<int>> &variant_position_dict,
                                                std::map<std::string, std::vector<std::string>> &variant_base_dict,
                                                std::map<std::string, std::vector<std::string>> &context_dict,
                                                std::map<std::string, int> &mapping_to_var,
                                                std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<int>>>> &matrix,
                                                std::map<int, std::vector<int>> &read_pos_loc,
                                                int HOMOPOLYMER_SIZE, const std::vector<std::string> &HOMOPOLYMER, int window);

void write_to_filter_cluster_bed(std::string filter_cluster_bed, std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<int>>>> &matrix,
                                 std::map<std::string, std::vector<int>> &variant_position_dict, std::map<int, std::vector<int>> &read_pos_loc,
                                 double CLUSTERED_FRAC, int NONCLUSTERED_COUNT);

#endif //INC_78_PY78MACHINELEARNEXTRACT_FEATURES_FROM_READS_FILTER_CLUSTER_H
