//
// Created by Colin Small on 8/2/21.
//

#ifndef CFSNV__JENKS_H
#define CFSNV__JENKS_H

#include <Rcpp.h>
using namespace Rcpp;

#include <tuple>
#include <vector>

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> jenksInitMatrices(std::vector<double> data, int nClasses);

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> jenksMatrices(std::vector<double> data, int nClasses);

std::vector<double> getJenksBreaks(std::vector<double> data, std::vector<std::vector<double>> lowerClassLimits, int nClasses);

std::vector<std::vector<double>> separateData(std::vector<double> data, std::vector<double> boundaries);

std::tuple<std::vector<std::vector<double>>, std::vector<double>> jenks(std::vector<double> data, int nClasses);

#endif //CFSNV__JENKS_H
