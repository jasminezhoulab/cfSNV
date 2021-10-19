//
// Created by Colin Small on 7/30/21.
//

#ifndef CFSNV_FILTER_WITH_PILEUP_H
#define CFSNV_FILTER_WITH_PILEUP_H

#include <Rcpp.h>
using namespace Rcpp;

#include <tuple>
#include <string>
#include <vector>

std::tuple<std::vector<std::string>, std::vector<bool>, std::string, bool, std::string, std::vector<bool>, std::string, bool, std::string, char, double> generateIntermediateResultOneLine(std::string line, double depth);

std::vector<double> generateIntermediateResultWholeFile(std::string filename, std::string output, std::string database, double depth);

std::vector<std::string> generateRecordOneLine(std::string line, double jenks_estimate);

std::tuple<std::vector<std::string>, std::string, std::string, std::string> decideOutputStatusOneLine(std::vector<std::string> record);

int generateRecordAndOutputWholeFile(std::string fintermediate, std::string foutputPass, std::string foutputCheck, std::string frecord, double jenksEstimate);

#endif //CFSNV_FILTER_WITH_PILEUP_H
