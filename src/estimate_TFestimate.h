//
// Created by Colin Small on 7/6/21.
//

#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "hotspot.h"

#ifndef CFSNV_ESTIMATE_TFESTIMATE_H
#define CFSNV_ESTIMATE_TFESTIMATE_H

void importPredefinedHotspotPrior(std::fstream file);

void selectHotSpot(std::string chrom,
                   std::string pos,
                   std::string basestring,
                   std::vector<double> quallist,
                   std::vector<double> maplist,
                   std::string basestringNormal,
                   std::vector<double> quallistNormal,
                   std::vector<double> maplistNormal,
                   std::string basestringExtendedFrags,
                   std::vector<double> quallistExtendedFrags,
                   std::vector<double> maplistExtendedFrags,
                   std::string basestringNotCombined,
                   std::vector<double> quallistNotCombined,
                   std::vector<double> maplistNotCombined,
                   std::map<char, int> basecountNotCombined,
                   std::map<char, int> basecountExtendedFrags,
                   char variantBase, int nread, double MERGED_VAF_THRESHOLD);

void selectHotspotFromPredefinedHotspotPrior(std::string chrom,
                                             std::string pos,
                                             std::string basestring,
                                             std::vector<double>& quallist,
                                             std::vector<double>& maplist,
                                             std::string basestringNormal,
                                             std::vector<double>& quallistNormal,
                                             std::vector<double>& maplistNormal,
                                             std::string basestringExtendedFrags,
                                             std::vector<double>& quallistExtendedFrags,
                                             std::vector<double>& maplistExtendedFrags,
                                             std::string basestringNotCombined,
                                             std::vector<double>& quallistNotCombined,
                                             std::vector<double>& maplistNotCombined,
                                             std::map<char, int>& basecountNotCombined,
                                             std::map<char, int>& basecountExtendedFrags,
                                             char variantBase);

double calculateTumorFractionLikelihood(double tumor_fraction,
                                     std::string basestring_merge,
                                     std::vector<double>& quallist_merge,
                                     std::vector<double>& maplist_merge,
                                     std::string basestring_normal,
                                     std::vector<double>& quallist_normal,
                                     std::vector<double>& maplist_normal,
                                     char variant_base);

std::tuple<int, std::vector<double>> estimateTumorFraction(hotspot hs);

#endif //CFSNV_ESTIMATE_TFESTIMATE_H
