//
// Created by Colin Small on 7/8/21.
//

#include <Rcpp.h>
using namespace Rcpp;

#include "hotspot.h"
#include <map>
#include <cmath>

std::vector<hotspot> HOTSPOT = {};
std::vector<std::tuple<std::string, int>> HOTSPOT_PRIOR = {};

std::map<std::string, double> PRIOR = {
{"AA/AA", std::pow(10.0,6.0)},
{"AA/AB", std::pow(10.0,4.0)},
{"AA/BB", std::pow(10.0,2.0)},
{"AB/AA", std::pow(10.0,2.0)},
{"AB/AB", std::pow(10.0,4.0)},
{"AB/BB", std::pow(10.0,2.0)},
{"BB/AA", 1.0},
{"BB/AB", 1.0},
{"BB/BB", std::pow(10.0,4.0)}
};

int GERMLINE_VARIANT_COUNT = 2;
double TRIALLELE_COUNT = 4;
double NORMAL_COUNT_BINOM = 0.05;
std::map<std::string, std::vector<std::string>> dbSNP;