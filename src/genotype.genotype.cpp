#include <Rcpp.h>
using namespace Rcpp;

#include "genotype.genotype.h"
#include "_parameter.h"
#include "_probability.h"
#include "_helper.h"
#include "_filter.h"

#include <float.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>

//# log( P( G, X | theta ) )
//# P( G | X , theta ) ~ P( G, X | theta )
std::map<std::string, double> calculateJointGenotypeLogposterior(double tumorFraction,
                                                                 std::string basestringMerge,
                                                                 std::vector<double> quallistMerge,
                                                                 std::vector<double> maplistMerge,
                                                                 std::string basestringNormal,
                                                                 std::vector<double> quallistNormal,
                                                                 std::vector<double> maplistNormal,
                                                                 char variantBase)
{
    //# input: tumor fraction, converted basestring, base quality list
    //# output: log posterior of all joint genotypes log( P( G, X | theta ) ), non-reference alleles, base count string
    //# P( G | X , theta ) is computed, but P( G | X , theta ) ~ P( G, X | theta )
    std::map<std::string, double> logposteriorDict = {};

    for(auto & it : PRIOR){
        std::string key = it.first;
        logposteriorDict[key] = (basestringMerge.size() + basestringNormal.size()) *
                ((1.0/basestringMerge.size()) *
                         calculateJointGenotypeTumorFractionLoglikelihood(tumorFraction, basestringMerge, quallistMerge,
                                                                          maplistMerge, key, variantBase) +
                  (1.0/basestringNormal.size()) *
                          calculateJointGenotypeTumorFractionLoglikelihood(0.0, basestringNormal, quallistNormal,
                                                                           maplistNormal, key, variantBase));
    }

    return logposteriorDict;
}

double calculateConfidence(std::map<std::string, double> logposteriorDict)
{
    //# input: log posterior list of all joint genotypes
    //# output: confidence score of the genotype with the maximum posterior
    //# confidence score is computed by the ratio of the maximum posterior and the second maximum posterior
    double maxVar = std::max(logposteriorDict["AA/AB"], logposteriorDict["AA/BB"]);
    std::vector<double> altList(0);
    for(auto& it : logposteriorDict){
        if(it.first != "AA/AB" and it.first != "AA/BB")
            altList.push_back(it.second);
    }
    double maxAlt = *std::max_element(altList.begin(), altList.end());
    double confidence = std::exp(maxVar)/std::exp(maxAlt);
    return confidence;
}

std::string calculateGenotype(std::map<std::string, double> logposteriorDict)
{
    //# input: log posterior list of all joint genotypes, non-reference alleles
    //# output: joint genotype, non-reference alleles
    //# joint genotype is determined by maximizing a posterior
    double maxVal = -DBL_MAX;
    std::string maxString;
    for(auto& it : logposteriorDict){
        if(it.second >= maxVal) {
            maxVal = it.second;
            maxString = it.first;
        }
    }
    return maxString;
}

void outputGenotyping(std::string outputFile,
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
                        std::vector<std::string> baseinfoList)
{
    //# input: chromosome, position, basestring, joint genotype, non-reference alleles, base count string, confidence, output file
    //# output: none
    //# only somatic variants are output currently
    //# all variants should be output
    //# related global variables:
    std::ofstream output;
    output.open(outputFile, std::ofstream::out | std::ofstream::app);

    output << "\"" << chrom << "\"" << "\t";
    output << "\"" << pos << "\"" << "\t";
    output << "\"" << (char)std::toupper(ref) << "\"" << "\t";
    output << "\"" << cntNormal << "\"" << "\t";
    output << "\"" << cntVarNormal << "\"" << "\t";
    output << "\"" << cntTumor << "\"" << "\t";
    output << "\"" << cntVarTumor << "\"" << "\t";
    output << "\"" << cntTumorMerge << "\"" << "\t";
    output << "\"" << cntVarTumorMerge << "\"" << "\t";
    output << "\"" << variantBase << "\"" << "\t";
    output << "\"" << jointGenotypeMAP << "\"" << "\t";
    output << "\"" << std::setprecision(12) << std::fixed << VAFTumor << "\"" << "\t";
    output << "\"" << std::setprecision(5) << std::fixed << confidence << "\"" << "\t";
    output << "\"" << jointGenotypeMAPMerge << "\"" << "\t";
    output << "\"" << std::setprecision(12) << std::fixed << VAFTumorMerge << "\"" << "\t";
    output << "\"" << std::setprecision(5) << std::fixed << confidenceMerge << "\"" << "\t";
    for(std::string& s : baseinfoList)
        output << "\"" << s << "\"" << "\t";
    output << "\n";

    output.close();
}

void callVariants(std::string finput,
                  std::string foutput_call,
                  double tumorFraction,
                  double depth,
                  double genotypeThreshold,
                  int GERMLINE_VARIANT_COUNT,
                  double NORMAL_COUNT_BINOM,
                  int MIN_PASS_SUPPORT_COUNT)
{
    std::ifstream inputFile(finput);

    // Clear existing output file text
    std::ofstream outputFile(foutput_call);
    outputFile.close();

    int n = 0;
    int m = 0;

    std::string lineText;
    while ( std::getline(inputFile, lineText)) {
        m += 1;
        std::vector<std::string> sp = split(lineText, "\t");
        sp[sp.size()-1] = trim_copy(sp[sp.size()-1]);
        std::string chrom = sp[0];
        int pos = stoi(sp[1]);
        char ref = sp[2][0];
        int nread = sp[19].size();
        std::string basestring = sp[19];
        std::string basestringAll = sp[4];

        if(basestring.empty())
            continue;

        std::vector<double> quallist = stringToQual(sp[20]);
        std::vector<double> maplist = stringToQual(sp[21]);

        std::string basestringNormal = sp[22];
        std::vector<double> quallistNormal = stringToQual(sp[23]);
        std::vector<double> maplistNormal = stringToQual(sp[24]);

        std::string baseStringExtendedFrags = sp[25];
        std::vector<double> quallistExtendedFrags = stringToQual(sp[26]);
        std::vector<double> maplistExtendedFrags = stringToQual(sp[27]);

        std::string basestringNotCombined = sp[28];
        std::vector<double> quallistNotCombined = stringToQual(sp[29]);
        std::vector<double> maplistNotCombined = stringToQual(sp[30]);

        std::string basestringNormalAll = sp[8];
        std::map<char, int> basecountNotCombinedAll = countBase(sp[12]);
        std::map<char, int> basecountExtendedFragsAll = countBase(sp[16]);
        char variantBaseAll = findMajorVariant(basecountNotCombinedAll, basecountExtendedFragsAll);

        std::map<char, int> basecountNotCombined = countBase(basestringNotCombined);
        std::map<char, int> basecountExtendedFrags = countBase(baseStringExtendedFrags);
        char variantBase = findMajorVariant(basecountNotCombined, basecountExtendedFrags);

        std::map<char, int> basecountTumorAll = countBase(basestringAll);
        std::map<char, int> basecountTumor = countBase(basestring);

        bool trialleleTumorAll = filterTriallelicPosition(basecountTumorAll, variantBaseAll, depth);
        bool trialleleTumor = filterTriallelicPosition(basecountTumor, variantBase, depth);

        if(variantBase != variantBaseAll)
            continue;

        if(std::count(basestring.begin(), basestring.end(), std::toupper(variantBase)) +
           std::count(basestring.begin(), basestring.end(), std::tolower(variantBase)) == 0)
            continue;

        if(!trialleleTumorAll or !trialleleTumor)
            continue;

        std::string basestringNormalAllUpper = toUpper(basestringNormalAll);
        if(std::count(basestringNormalAllUpper.begin(), basestringNormalAllUpper.end(), variantBaseAll) > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER)
            continue;
        if(std::count(basestringNormalAllUpper.begin(), basestringNormalAllUpper.end(), variantBase) > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER)
            continue;

        std::string basestringNormalUpper = toUpper(basestringNormal);
        if(std::count(basestringNormalUpper.begin(), basestringNormalUpper.end(), variantBase) > GERMLINE_VARIANT_COUNT)
            continue;

        std::string basestringMerge = baseStringExtendedFrags + basestringNotCombined;

        std::vector<double> quallistMerge;
        quallistMerge.reserve(quallistExtendedFrags.size() + quallistNotCombined.size());
        quallistMerge.insert(quallistMerge.end(), quallistExtendedFrags.begin(), quallistExtendedFrags.end());
        quallistMerge.insert(quallistMerge.end(), quallistNotCombined.begin(), quallistNotCombined.end());

        std::vector<double> maplistMerge;
        maplistMerge.reserve(maplistExtendedFrags.size() + maplistNotCombined.size());
        maplistMerge.insert(maplistMerge.end(), maplistExtendedFrags.begin(), maplistExtendedFrags.end());
        maplistMerge.insert(maplistMerge.end(), maplistNotCombined.begin(), maplistNotCombined.end());

        if(basestringNormal.size() < DEPTH_FOR_DETECTION_NORMAL)
            continue;
        if(basestringMerge.size() < DEPTH_FOR_DETECTION_TUMOR or basestring.size() < DEPTH_FOR_DETECTION_TUMOR)
            continue;

        std::string basestringMergeUpper = toUpper(basestringMerge);
        std::string basestringUpper = toUpper(basestring);
        double cntVarTumor = std::count(basestringUpper.begin(), basestringUpper.end(), variantBase);
        double cntVarTumorMerge = std::count(basestringMergeUpper.begin(), basestringMergeUpper.end(), variantBase);
        double cntVarNormal = std::count(basestringNormalUpper.begin(), basestringNormalUpper.end(), variantBase);

        double VAFTumorMerge = cntVarTumorMerge/basestringMerge.size();
        double VAFTumor = cntVarTumor/basestring.size();
        double VAFNormal = cntVarNormal/basestringNormal.size();

        int cntAltNormal = basestringNormalUpper.size() - std::count(basestringNormalUpper.begin(), basestringNormalUpper.end(), 'R');
        int cntNormal = basestringNormal.size();

        int NORMAL_COUNT_VAR = floor(0.5*MIN_PASS_SUPPORT_COUNT);
        if(VAFNormal > SOMATIC_VAF_THRESHOLD_IN_NORMAL)
            continue;
        if(cntVarNormal > 0 and binomialPMF(cntVarNormal, basestringNormal.size(), VAFTumor) > NORMAL_COUNT_BINOM)
            continue;
        if(cntVarNormal > NORMAL_COUNT_VAR or cntAltNormal > NORMAL_COUNT_ALT)
            continue;
        if(cntVarNormal > 0 and (VAFNormal >= VAFTumorMerge or VAFNormal >= VAFTumor or cntVarNormal >= cntVarTumor or cntVarNormal >= cntVarTumorMerge))
            continue;

        std::map<std::string, double> logposteriorDictMerge = calculateJointGenotypeLogposterior(tumorFraction, basestringMerge, quallistMerge, maplistMerge, basestringNormal, quallistNormal, maplistNormal, variantBase);
        double confidenceMerge = calculateConfidence(logposteriorDictMerge);
        std::string jointGenotypeMAPMerge = calculateGenotype(logposteriorDictMerge);

        std::map<std::string, double> logposteriorDict = calculateJointGenotypeLogposterior(tumorFraction, basestring, quallist, maplist, basestringNormal, quallistNormal, maplistNormal, variantBase);
        double confidence = calculateConfidence(logposteriorDict);
        std::string jointGenotypeMAP = calculateGenotype(logposteriorDict);

        std::vector<std::string>subsp(sp.begin() + 19, sp.begin() + 31);


        if(jointGenotypeMAP == "AA/BB" or jointGenotypeMAP == "AA/AB" or jointGenotypeMAPMerge == "AA/BB" or jointGenotypeMAPMerge == "AA/AB"){
            n += 1;
            if( VAFTumor <= genotypeThreshold and VAFTumorMerge <= genotypeThreshold)
                outputGenotyping(foutput_call, chrom, pos, ref, cntNormal, cntVarNormal, basestring.size(), cntVarTumor, basestringMerge.size(), cntVarTumorMerge, variantBase, jointGenotypeMAP, VAFTumor, confidence, jointGenotypeMAPMerge, VAFTumorMerge, confidenceMerge, subsp);
        }
    }
    inputFile.close();
}

// [[Rcpp::export]]
int genotype_genotype_main(std::string finput, std::string foutputCall, double tumorFraction,
                           double mergedVAFThreshold, double depth, int MIN_PASS_SUPPORT_COUNT){
    // std::string finput = argv[1];
    // std::string foutputCall = argv[2];
    // double tumorFraction = std::stod(argv[3]);
    // double mergedVAFThreshold = std::stod(argv[4]);
    // double depth = std::stod(argv[5]);

    double genotypeThreshold = mergedVAFThreshold * 1.2;
    GERMLINE_VARIANT_COUNT = getGermlineVariantCountThreshold(depth);
    NORMAL_COUNT_BINOM = getNormalCountBinomThreshold(depth);
    callVariants(finput, foutputCall, tumorFraction, depth, genotypeThreshold, GERMLINE_VARIANT_COUNT, NORMAL_COUNT_BINOM, MIN_PASS_SUPPORT_COUNT);
    return 0;
}
