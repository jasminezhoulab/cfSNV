//
// Created by Colin Small on 7/7/21.
//

#ifndef CFSNV__FILTER_H
#define CFSNV__FILTER_H

#include <Rcpp.h>
using namespace Rcpp;

#include "string"
#include "map"
#include <vector>

bool filterStrandBiasMerge(std::map<char, int> & basecountNotCombined,
                           std::map<char, int> & basecountExtendedFrags,
                           char variantBase);

bool filterStrandBiasUnmerge(std::map<char, int> & basecount,
                             char variantBase);

std::tuple<bool, bool> filterBothStrandAboveAverageMerge(std::map<char, int> & basecountNotCombined,
                                       std::map<char, int> & basecountExtendedFrags);


bool filterSupport(std::string basestringMerge, char variantBase);

std::tuple<bool, bool> filterBothStrandAboveAverageUnmerge(std::map<char, int> basecount);

int getGermlineVariantCountThreshold(double depth);

char findMajorVariant(std::map<char, int> & basecountNC, std::map<char, int> & basecountOF);

bool filterTriallelicPosition(std::map<char, int> basecount, char var, int depth);

double getNormalCountBinomThreshold(double depth);

bool filterAllMappingQuality(std::vector<double> maplist);

bool filterVariantMappingQuality(std::vector<double> mapList, std::string basestring, char variantBase);

bool filterVariantBaseQuality(std::vector<double> quallist, std::string basestring, char variantBase);

bool filterSupportingReadsUnmerge(std::string basestring, char variantBase);

std::tuple<bool, std::string> filterVAF(double vaf, double depth);

bool filterNormalCoverage(int countNormal);

std::tuple<bool, std::string> filterTumorCoverage(int count);

void importdbSNP(std::string dbSPNFile);

bool filterSupportingCountByBinomialTest(double jenksEstimate, std::string basestring, char variantBase);

std::tuple<double, int, int> finalEstimationWithJenks(std::vector<double> VAFList);

#endif //CFSNV__FILTER_H
