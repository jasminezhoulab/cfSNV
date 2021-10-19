//
// Created by Colin Small on 7/6/21.
//

#ifndef CFSNV__PARAMETER_H
#define CFSNV__PARAMETER_H

#include <Rcpp.h>
using namespace Rcpp;

#include "hotspot.h"
#include <vector>
#include <map>
#include <cmath>

const double BASEQUAL_THRESHOLD = 0.05;
const double MAPQUAL_THRESHOLD = 0.3;
const int COUNT_VAR_TUMOR_HIGHQUAL = 3;


// py2
const int DEPTH_FOR_GERMLINE = 30;
const double GERMLINE_VAF_THRESHOLD_IN_NORMAL = 0.01;
const int GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER = 2;
extern int GERMLINE_VARIANT_COUNT;
// const int MAX_NUMBER_OF_HOTPOST = 50;
const double AVG_MAPQUAL_FOR_ESTIMATION = 0.1;
const int ZERO_MAPQUAL_COUNT_FOR_ESTIMATION = 5;

const double GRIDWIDTH = 0.001;
const int MAXSEARCH = 400;
const int MAXSEARCH_WITH_NORMAL = 1000;
const int DEPTH_FOR_ESTIMATION = 80;

// _filter
const double STRAND_BIAS_RATIO_VARIANT_TO_ALL = 5.0;
const double STRAND_BIAS_BINOMIAL_PROB = 0.05;
const double THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF = 0.8;
const int COUNT_FOR_STRONG_EVIDENCE = 3;

const double TRIALLELE_VAF_RATIO = 0.5;
const double TRIALLELE_VAF = 0.02;
extern double TRIALLELE_COUNT;

extern std::vector<hotspot> HOTSPOT;
extern std::vector<std::tuple<std::string, int>> HOTSPOT_PRIOR;

const std::map<std::string, double> EST_PRIOR = {
        {"AA/AA", std::pow(10.0,6.0)},
        {"AA/AB", std::pow(10.0,4.0)},
        {"AA/BB", std::pow(10.0,2.0)},
        {"AB/AA", std::pow(10.0,2.0)},
        {"AB/AB", std::pow(10.0,3.0)},
        {"AB/BB", 0.0001},
        {"BB/AA", 0.0001},
        {"BB/AB", 0.0001},
        {"BB/BB", 0.0001}
};

// py3

const int NORMAL_COUNT_ALT = 7;
//const int NORMAL_COUNT_VAR = 2;
const double SOMATIC_VAF_THRESHOLD_IN_NORMAL = 0.12;
extern double NORMAL_COUNT_BINOM;
const int NORMAL_COUNT_FOR_VAF_THRESHOLD = 50;
const int DEPTH_FOR_DETECTION_NORMAL = 4;
const int DEPTH_FOR_DETECTION_TUMOR = 4;

extern std::map<std::string, double> PRIOR;

//

const double VAF_FOR_STRONG_EVIDENCE = 0.0005;
const int LOW_MAPQUAL_COUNT_FOR_DETECTION_VARIANT = 3;
const int DEPTH_ACTIVATE_LOW_MAPQUAL_COUNT = 6;
const int LOW_MAPQUAL_COUNT_FOR_DETECTION_ALL = 10;
const int DEPTH_ACTIVATE_LOW_MAPQUAL_COUNT_ALL = 20;
const double LOW_MAPQUAL = 0.3;
const double OK_MAPQUAL = 0.1;
const double OK_BASEQUAL = 0.005;
const int THRESHOLD_OK_BASEQUAL_COUNT = 3;
const int THRESHOLD_OK_MAPQUAL_COUNT = 3;
const int THRESHOLD_OK_MAPQUAL_COUNT_ALL = 20;
const double AVG_MAPQUAL_FOR_DETECTION = 0.1;
const double LOW_MAPQUAL_FRAC_FOR_DETECTION = 0.4;
const double THRESHOLD_VARIANT_QUALITY_TO_REFERENCE_QUALITY = 7.0;
const int VARIANT_COUNT_FOR_JENK_ESTIMATION = 30;
const double P_VALUE_BINOMIAL_TEST_JENKS_ESTIMATE = 0.1;
const int VARIANT_COUNT_STOP_ADDING_IN_JENK_ESTIMATION = 20;
const std::string CIGAR_CHARACTER = "MSIDHNPX=";
const int DISTANCE_FROM_INDEL = 4;
const int NUMBER_INDEL_NEARBY = 3;
const int THRESHOLD_ADJACENT_ALL_T_STATISTICS = 10;
const int THRESHOLD_ADJACENT_VARIANT_T_STATISTICS = 10;

const int MIN_HOLD_SUPPORT_COUNT = 12;
const int MIN_PASS_SUPPORT_COUNT = 5;

extern std::map<std::string, std::vector<std::string>> dbSNP;

#endif //CFSNV__PARAMETER_H
