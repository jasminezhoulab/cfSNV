//
// Created by Colin Small on 7/19/21.
//

#include <Rcpp.h>
using namespace Rcpp;

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>

#ifndef CFSNV_GENOTYPE_GENOTYPE_H
#define CFSNV_GENOTYPE_GENOTYPE_H

std::map<std::string, double> calculateJointGenotypeLogposterior(double tumorFraction,
                                                                 std::string basestringMerge,
                                                                 std::vector<double> quallistMerge,
                                                                 std::vector<double> maplistMerge,
                                                                 std::string basestringNormal,
                                                                 std::vector<double> quallistNormal,
                                                                 std::vector<double> maplistNormal,
                                                                 char variantBase);

double calculateConfidence(std::map<std::string, double> logposteriorDict);

std::string calculateGenotype(std::map<std::string, double> logposteriorDict);

void outputGenotyping(std::string output,
                        std::string chrom,
                        int pos,
                        char ref,
                        int cntNormal,
                        int cntVarNormal,
                        int cntTumor,
                        int cntVarTumor,
                        int cntTumorMerge,
                        int cntVarTumorMerge,
                        char variantBase,
                        std::string jointGenotypeMAP,
                        double VAFTumor,
                        double confidence,
                        std::string jointGenotypeMAPMerge,
                        double VAFTumorMerge,
                        double confidenceMerge,
                        std::vector<std::string> baseinfoList);

void callVariants(std::string finput,
                  std::string foutput_call,
                  double tumorFraction,
                  int depth,
                  double genotypeThreshold,
                  int GERMLINE_VARIANT_COUNT,
                  double NORMAL_COUNT_BINOM,
                  int MIN_PASS_SUPPORT_COUNT);


#endif //CFSNV_GENOTYPE_GENOTYPE_H
