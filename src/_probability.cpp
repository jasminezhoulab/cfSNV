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
#include <boost/math/distributions/binomial.hpp>


double binomialPMF(int k, int n, double p){
//    double binomialCoefficient = (double) factorial(n) / (factorial(k) * factorial(n-k));
//    if( isnan(binomialCoefficient))
//        throw std::overflow_error("Factorial overflow in binomialPMF()");
//    return binomialCoefficient * std::pow(p, k) * std::pow(1-p, n-k);
    //std::cout << n << " " << p << std::endl;
    return boost::math::pdf(boost::math::binomial(n, p), k);
}

//double binomialCDF(int k, int n, double p){
//    double cdf = 0;
//
//    for(int i = 0; i < k+1; i++){
//        double binomialCoefficient = factorial(n) / (factorial(i) * factorial(n-i));
//        cdf += binomialCoefficient * std::pow(p, i) * std::pow(1-p, n-i);
//    }
//    return (double) cdf;
//}


double binomialCDF(int k, int n, double p){
    return boost::math::cdf(boost::math::binomial(n, p), k);
}



std::map<char, int> countBase(std::string basestring) {
    std::map<char, int> basecount = {{'A',0}, {'a',0}, {'C',0}, {'c',0}, {'G',0}, {'g',0}, {'T',0}, {'t',0}, {'R',0}, {'r',0}};

    for(char c: basestring){
        basecount[c] += 1;
    }

    return basecount;
}

void translate_basestring(std::string & basestring, char & variant_base) {
    // Variant_base must be lowercase
    std::replace( basestring.begin(), basestring.end(), variant_base, 'X');
    std::replace( basestring.begin(), basestring.end(), (char)std::tolower(variant_base), 'X');
    std::replace( basestring.begin(), basestring.end(), 'r', 'R');
}

std::map<char, double> observeVariantProbability(double tumorFraction, std::string joint_genotype) {
    std::map<char, double> variantAlleleObservedProbability = {};

    // TODO: Consider hashing the string and using a switch statement

    if(joint_genotype == "AA/AA") {
        variantAlleleObservedProbability['R'] = 1.0;
        variantAlleleObservedProbability['X'] = 0.0;
    }
    else if(joint_genotype == "AA/AB") {
        variantAlleleObservedProbability['R'] = 1.0 - 0.5 * tumorFraction;
        variantAlleleObservedProbability['X'] = 0.5 * tumorFraction;
    }
    else if(joint_genotype == "AA/BB") {
        variantAlleleObservedProbability['R'] = 1.0 - tumorFraction;
        variantAlleleObservedProbability['X'] = tumorFraction;
    }
    else if(joint_genotype == "AB/AA") {
        variantAlleleObservedProbability['R'] = 0.5 + 0.5 * tumorFraction;
        variantAlleleObservedProbability['X'] = 0.5 - 0.5 * tumorFraction;
    }
    else if(joint_genotype == "AB/AB"){
        variantAlleleObservedProbability['R'] = 0.5;
        variantAlleleObservedProbability['X'] = 0.5;
    }
    else if(joint_genotype == "AB/BB") {
        variantAlleleObservedProbability['R'] = 0.5 - 0.5 * tumorFraction;
        variantAlleleObservedProbability['X'] = 0.5 + 0.5 * tumorFraction;
    }
    else if(joint_genotype == "BB/AA") {
        variantAlleleObservedProbability['R'] = tumorFraction;
        variantAlleleObservedProbability['X'] = 1.0 - tumorFraction;
    }
    else if(joint_genotype == "BB/AB") {
        variantAlleleObservedProbability['R'] = 0.5 * tumorFraction;
        variantAlleleObservedProbability['X'] = 1.0 - 0.5 * tumorFraction;
    }
    else if(joint_genotype == "BB/BB") {
        variantAlleleObservedProbability['R'] = 0.0;
        variantAlleleObservedProbability['X'] = 1.0;
    }
    else
        std::cerr << "Error with input genotype";

    return variantAlleleObservedProbability;
}

double calculateJointGenotypeTumorFractionLoglikelihood(double tumor_fraction,
                                                        std::string& basestring,
                                                        std::vector<double>& quallist,
                                                        std::vector<double>& maplist,
                                                        std::string& joint_genotype,
                                                        char& variant_base) {
    translate_basestring(basestring, variant_base);
    std::vector<int> rBasestringInterestIndices = getIndicesWhereEqualFromString(basestring, 'R');
    std::vector<int> xBasestringInterestIndices = getIndicesWhereEqualFromString(basestring, 'X');

    rBasestringInterestIndices.insert(rBasestringInterestIndices.end(), xBasestringInterestIndices.begin(), xBasestringInterestIndices.end());
    std::vector<int> interestIndices = rBasestringInterestIndices;
    std::sort(interestIndices.begin(), interestIndices.end());
    //    std::vector<int> interestIndices = rBasestringInterestIndices + xBasestringInterestIndices;
    std::vector<double> interestQuallist = accessMultipleIndicesVector(quallist, interestIndices);
    std::vector<double> interestMaplist = accessMultipleIndicesVector(maplist, interestIndices);

    std::map<char, double> variantAlleleFrequences = observeVariantProbability(tumor_fraction, joint_genotype);
    std::vector<char> interestBasestring = accessMultipleIndicesString(basestring, interestIndices);

    std::vector<double> interestVAFList;
    interestVAFList.reserve(interestBasestring.size());
    for(char c : interestBasestring){
        interestVAFList.push_back(variantAlleleFrequences[c]);
    }

    std::vector<double> onevec(interestQuallist.size(), 1.0);
    std::vector<double> loglist = veclog(interestVAFList*(onevec-interestQuallist)*(onevec-interestMaplist) + (onevec-interestVAFList)*(interestQuallist+interestMaplist+interestQuallist*interestMaplist));

    double loglikelihood = 0.0;
    for(double& d : loglist){
        loglikelihood += d;
    }

    return loglikelihood;
}

std::vector<double> stringToQual(std::string qualityString){
    std::vector<double> qualityVector;
    qualityVector.reserve(qualityString.size());

    for(char& c : qualityString){
        double y = (0-((double)c-33.0)/10.0);
        qualityVector.push_back(std::pow(10.0, y));
    }

    return qualityVector;
}
