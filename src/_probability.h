//
// Created by Colin Small on 7/6/21.
//

#ifndef CFSNV__PROBABILITY_H
#define CFSNV__PROBABILITY_H

#include <Rcpp.h>
using namespace Rcpp;

#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <valarray>
//#include <string_view>
#include <iostream>
#include <chrono>
#include "_helper.h"


std::map<char, int> countBase(std::string basestring);


/**
 * Convert reference alleles to 'R' and non-reference alleles to 'X'
 *
 * @param basestring Basestring in rebase file
 * @param variant_base
 * @return
 */
void translate_basestring(std::string& basestring, char& variant_base);



std::map<char, double> observeVariantProbability(double tumorFraction, std::string joint_genotype);



double calculateJointGenotypeTumorFractionLoglikelihood(double tumor_fraction,
                                                        std::string& basestring,
                                                        std::vector<double>& quallist,
                                                        std::vector<double>& maplist,
                                                        std::string& joint_genotype,
                                                        char& variant_base
);

double binomialPMF(int k, int n, double p);

double binomialCDF(int k, int n, double p);

std::vector<double> stringToQual(std::string qualityString);


//int main(){
//    auto start = std::chrono::high_resolution_clock::now();
//    std::string basestring = "rrRrrRRRrRRrrRRRRrrrrRrRrrRrRRrRaaaAAAAAAaaAaAaaaAaARrRRRrrrRrRrrRrRRRRrRrAAaaaaAaaaAAaaAaa";
//    std::string joint_genome = "BB/BB";
//    std::vector<double> maplist = {0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.251188644, 0.000001, 0.000001, 0.000001};
//    std::vector<double> quallist = {0.022535235, 0.00528222, 0.009541663, 0.013264347, 0.011090627, 0.003555052, 0.027080566, 0.015750399, 0.022683902, 0.031131367, 0.018204997, 0.00997583, 0.005136876, 0.019508333, 0.028207442, 0.028061324, 0.027953896, 0.028271181, 0.013522274, 0.026801979, 0.0016488, 0.005067841, 0.031046452, 0.019417774, 0.027828293, 0.023785416, 0.019658524, 0.006075595, 0.007924912, 0.029089328, 0.010350578, 0.017759347, 0.007546654, 0.021277069, 0.018700052, 0.006809801, 0.027671635, 0.006925604, 0.028333594, 0.024216763, 0.028855102, 0.029491, 0.011219381, 0.022957093, 0.014945312, 0.011761302, 0.001022148, 0.019425203, 0.026984942, 0.005588313, 0.028931969, 0.014299779, 0.00322236, 0.021123591, 0.004383596, 0.02223039, 0.025926345, 0.020896607, 0.019617726, 0.005057836, 0.012143053, 0.026517415, 0.020205884, 0.00856174, 0.007850494, 0.019941156, 0.028486627, 0.020179196, 0.013996312, 0.026127624, 0.01342404, 0.012679071, 0.011769746, 0.021126596, 0.013215724, 0.00816966, 0.018582238, 0.004288396, 0.014510591, 0.022444977, 0.020117276, 0.005726826, 0.003265246, 0.023436632, 0.024888499, 0.028886556, 0.026453318, 0.023166153, 0.017933538, 0.024731459, 0.03156874};
//    double tumor_fraction = 0.0;
//    char variant_base = 'A';
//
//    for(int i = 0; i < 10000; i++) {
//        calculate_joint_genotype_tumor_fraction_loglikelihood(tumor_fraction,
//                                                              basestring,
//                                                              quallist,
//                                                              maplist,
//                                                              joint_genome,
//                                                              variant_base);
//    }
//    auto end = std::chrono::high_resolution_clock::now();
//    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()/1000000.0;
//}

#endif //CFSNV__PROBABILITY_H
