//
// Created by Colin Small on 7/7/21.
//
#include <Rcpp.h>
using namespace Rcpp;

#include "_filter.h"
#include "_parameter.h"
#include "_probability.h"
#include "_helper.h"
#include "_jenks.h"

#include <random>
#include <iostream>
#include <fstream>

bool filterStrandBiasMerge(std::map<char, int> & basecountNotCombined,
                                  std::map<char, int> & basecountExtendedFrags,
                                  char variantBase){
    int forward = basecountNotCombined['R'] + basecountNotCombined['A'] + basecountNotCombined['C'] + basecountNotCombined['T'] + basecountNotCombined['G'] + (basecountExtendedFrags['r'] + basecountExtendedFrags['a'] + basecountExtendedFrags['c'] + basecountExtendedFrags['t'] + basecountExtendedFrags['g'] + basecountExtendedFrags['R'] + basecountExtendedFrags['A'] + basecountExtendedFrags['C'] + basecountExtendedFrags['T'] + basecountExtendedFrags['G']);
    int reverse = basecountNotCombined['r'] + basecountNotCombined['a'] + basecountNotCombined['c'] + basecountNotCombined['t'] + basecountNotCombined['g'] + (basecountExtendedFrags['R'] + basecountExtendedFrags['A'] + basecountExtendedFrags['C'] + basecountExtendedFrags['T'] + basecountExtendedFrags['G'] + basecountExtendedFrags['r'] + basecountExtendedFrags['a'] + basecountExtendedFrags['c'] + basecountExtendedFrags['t'] + basecountExtendedFrags['g']);
    int forwardVar = basecountNotCombined[variantBase] + basecountExtendedFrags[std::tolower(variantBase)] + basecountExtendedFrags[variantBase];
    int reverseVar = basecountNotCombined[std::tolower(variantBase)] + basecountExtendedFrags[variantBase] + basecountExtendedFrags[std::tolower(variantBase)];

    if(reverseVar == 0 && forwardVar == 0)
        return false;
    if(reverseVar == 0 && reverse == 0)
        return true;
    if(reverseVar == 0 and reverse != 0) {
        if ( binomialPMF(forwardVar, forwardVar + reverseVar, (double) forward / (forward + reverse)) <
             STRAND_BIAS_BINOMIAL_PROB )
            return false;
        else
            return true;
    }
    else if(forwardVar == 0 and forward == 0)
        return true;
    else if(forwardVar == 0 and forward != 0)
    {
        if ( binomialPMF(reverseVar, forwardVar + reverseVar, (double) reverse / (forward + reverse)) <
             STRAND_BIAS_BINOMIAL_PROB )
            return false;
        else
            return true;
    }
    else{
        double varRatio = (double)forwardVar/reverseVar;
        double ratio = (double)forward/reverse;
        if(varRatio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL ||
           ratio/varRatio > STRAND_BIAS_RATIO_VARIANT_TO_ALL)
            return false;
        else
            return true;
    }
}

std::tuple<bool, bool> filterBothStrandAboveAverageMerge(std::map<char, int> & basecountNotCombined,
                                       std::map<char, int> & basecountExtendedFrags)
{
    char bothObserved = false;
    char aboveAverage = false;

    for( char base : {'A', 'C', 'T', 'G'})
    {
        if((basecountNotCombined[std::tolower(base)] > 0 && basecountNotCombined[base] > 0 ) ||
           basecountExtendedFrags[std::tolower(base)] > 0 ||
           basecountExtendedFrags[base] > 0)
            bothObserved = true;

        if( basecountNotCombined[base] + basecountNotCombined[std::tolower(base)] + 2*basecountExtendedFrags[std::tolower(base)] + 2*basecountExtendedFrags[base] >=
            THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF * (sumMap(basecountNotCombined) + 2*sumMap(basecountExtendedFrags) - basecountNotCombined['R'] - basecountNotCombined['r'] - 2*basecountExtendedFrags['R'] - 2*basecountExtendedFrags['r']))
            aboveAverage = true;
    }

    return std::tuple<bool, bool>{bothObserved, aboveAverage};
}

bool filterSupport(std::string basestring, char variantBase){
    int variantCount = std::count(basestring.begin(), basestring.end(), variantBase);
    int variantCountLower = std::count(basestring.begin(), basestring.end(), std::tolower(variantBase));
    if( variantCount + variantCountLower < COUNT_FOR_STRONG_EVIDENCE )
        return false;
    else
        return true;
}

bool filterStrandBiasUnmerge(std::map<char, int> & basecount,
                             char variantBase)
{
    int forward = basecount['R'] + basecount['A'] + basecount['C'] + basecount['T'] + basecount['G'];
    int reverse = basecount['r'] + basecount['a'] + basecount['c'] + basecount['t'] + basecount['g'];
    int forwardVar = basecount[variantBase];
    int reverseVar = basecount[std::tolower(variantBase)];

    if(reverseVar == 0 and reverse == 0)
        return true;
    if(reverseVar == 0 and reverse != 0)
    {
        if(binomialPMF(forwardVar, forwardVar+reverseVar, (double)forward/(forward+reverse)) < STRAND_BIAS_BINOMIAL_PROB)
            return false;
        else
            return true;
    }
    if(forwardVar == 0 and forward == 0)
        return true;
    if(forwardVar == 0 and forward != 0)
    {
        if( binomialPMF(reverseVar, forwardVar+reverseVar, (double)reverse/(forward+reverse)) < STRAND_BIAS_BINOMIAL_PROB)
            return false;
        else
            return true;
    }
    else{
        double varRatio = (double)forwardVar/reverseVar;
        double ratio = (double)forward/reverse;
        if(varRatio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL ||
           ratio/varRatio > STRAND_BIAS_RATIO_VARIANT_TO_ALL)
            return false;
        else
            return true;
    }
}

std::tuple<bool, bool> filterBothStrandAboveAverageUnmerge(std::map<char, int> basecount){
    char bothObserved = false;
    char aboveAverage = false;

    for(char base : {'A', 'C', 'T', 'G'})
    {
        if(basecount[std::tolower(base)] > 0 and basecount[base] > 0)
            bothObserved = true;
        if(basecount[base] + basecount[std::tolower(base)] >= THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF * (sumMap(basecount) - basecount['R'] - basecount['r']))
            aboveAverage = true;
    }

    return std::tuple<bool, bool>{bothObserved, aboveAverage};
}

int getGermlineVariantCountThreshold(double depth){
    return (int)std::log10(depth)-1;
}

char findMajorVariant(std::map<char, int> & basecountNC, std::map<char, int> & basecountOF){
    std::vector<char> alternativeNucleotide{'A', 'C', 'T', 'G'};
    std::vector<int> count;
    count.reserve(alternativeNucleotide.size());
    for(char& nucleotide : alternativeNucleotide){
        count.push_back(basecountNC[nucleotide] + basecountNC[std::tolower(nucleotide)] + 2*basecountOF[nucleotide] + 2*basecountOF[std::tolower(nucleotide)]);
    }
    int variantBaseIndex = std::max_element(count.begin(), count.end()) - count.begin();
    return alternativeNucleotide[variantBaseIndex];
}

bool filterTriallelicPosition(std::map<char, int> basecount, char var, int depth){
    int varCount = basecount[var] + basecount[std::tolower(var)];
    int totalCount = sumMap(basecount);
//    if(totalCount == 0)
//        std::cout << totalCount << "," << varCount << "," << basecount << ',' << var;
    double varVAF = (double)varCount/totalCount;
    TRIALLELE_COUNT = (int)std::log10(depth) + 2;
    int otherCount = 0;

    for(char c : "ACGT"){
        if(c == var)
            continue;
        int tmpCount = basecount[c] + basecount[std::tolower(c)];
        if( tmpCount > otherCount)
            otherCount = tmpCount;
    }

    double otherVAF = (double)otherCount/totalCount;
    if(varCount == 0 or otherVAF/varVAF > TRIALLELE_VAF_RATIO)
        return false;
    if(otherVAF > TRIALLELE_VAF and totalCount >= 80)
        return false;
    if(otherCount > TRIALLELE_COUNT and totalCount <= 140)
        return false;
    return true;
}

double getNormalCountBinomThreshold(double depth){
    if(depth < 1000)
        return 0.05;
    else
        return 0.2;
}

bool filterAllMappingQuality(std::vector<double> maplist){
    std::vector<int> mapListLessThanOKMAPQUALIndices = getIndicesWhereLessThan(maplist, OK_MAPQUAL);
    std::vector<int> maplistGreaterThanLOWMAPQUALIndices = getIndicesWhereGreaterThan(maplist, LOW_MAPQUAL);

    if(mapListLessThanOKMAPQUALIndices.size() > THRESHOLD_OK_MAPQUAL_COUNT_ALL)
        return true;

    if((maplistGreaterThanLOWMAPQUALIndices.size() > LOW_MAPQUAL_COUNT_FOR_DETECTION_ALL and maplist.size() < DEPTH_ACTIVATE_LOW_MAPQUAL_COUNT_ALL)
        or (accessMultipleIndicesVector(maplist, maplistGreaterThanLOWMAPQUALIndices).size()/maplist.size() > LOW_MAPQUAL_FRAC_FOR_DETECTION))
        return false;

    if( vecMedian(maplist) > AVG_MAPQUAL_FOR_DETECTION)
        return false;

    else
        return true;
}

bool filterVariantMappingQuality(std::vector<double> mapList, std::string basestring, char variantBase){
    std::string baseStringUpper = toUpper(basestring);
    std::vector<int> variantMaplistIndices = getIndicesWhereEqualFromString(baseStringUpper, variantBase);
    std::vector<double> variantMapList = accessMultipleIndicesVector(mapList, variantMaplistIndices);

    std::vector<double> variantMapListLessThanOKMAPQUAL = accessMultipleIndicesVector(variantMapList, getIndicesWhereLessThanOrEqual(variantMapList, OK_MAPQUAL));
    std::vector<double> variantMapListGreaterThanLOWMAPQUAL = accessMultipleIndicesVector(variantMapList,
                                                                                          getIndicesWhereGreaterThan(variantMapList, LOW_MAPQUAL));

//    double variantMaplistSum = sumVector(variantMapList);

    if(variantMapList.size() == 0)
        return false;
    if( variantMapListLessThanOKMAPQUAL.size() > THRESHOLD_OK_MAPQUAL_COUNT)
        return true;
    if( variantMapListGreaterThanLOWMAPQUAL.size() > LOW_MAPQUAL_COUNT_FOR_DETECTION_VARIANT and variantMapList.size() < DEPTH_ACTIVATE_LOW_MAPQUAL_COUNT)
        return false;
    if( variantMapListGreaterThanLOWMAPQUAL.size()/variantMapList.size() > LOW_MAPQUAL_FRAC_FOR_DETECTION)
        return false;
    if( vecMedian(variantMapList) > AVG_MAPQUAL_FOR_DETECTION)
        return false;
    return true;
}

bool filterVariantBaseQuality(std::vector<double> quallist, std::string basestring, char variantBase){
    std::string baseStringUpper = toUpper(basestring);
    std::vector<int> varID = getIndicesWhereEqualFromString(baseStringUpper, variantBase);
    std::vector<int> refID = getIndicesWhereEqualFromString(baseStringUpper, 'R');

    if(varID.size() == 0)
        return false;
    if(refID.size() == 0)
        return true;

    double avgVarQual = vecMedian(accessMultipleIndicesVector(quallist, varID));
    double avgRefQual = vecMedian(accessMultipleIndicesVector(quallist, refID));

    std::vector<double> quallistVarID = accessMultipleIndicesVector(quallist, varID);
    std::vector<int> quallistLessThanOKBASEQUAL = getIndicesWhereLessThan(quallistVarID, OK_BASEQUAL);

    if(avgVarQual/avgRefQual > THRESHOLD_VARIANT_QUALITY_TO_REFERENCE_QUALITY and quallistLessThanOKBASEQUAL.size() < THRESHOLD_OK_BASEQUAL_COUNT)
        return false;
    else
        return true;
}

bool filterSupportingReadsUnmerge(std::string basestring, char variantBase){
    int basecount = std::count(basestring.begin(), basestring.end(), variantBase);
    int basecountUpper = std::count(basestring.begin(), basestring.end(), std::tolower(variantBase));

    if(basecount + basecountUpper < COUNT_FOR_STRONG_EVIDENCE)
        return false;
    else
        return true;
}

std::tuple<bool, std::string> filterVAF(double vaf, double depth){
    double VAFForStrongEvidence = 3.0/depth;
    if(vaf < VAFForStrongEvidence)
        return std::tuple<char, std::string>(false, "weak VAF");
    else
        return std::tuple<char, std::string>(true, "strong VAF");
}

bool filterNormalCoverage(int countNormal){
    if(countNormal <= 7)
        return false;
    else
        return true;
}

std::tuple<bool, std::string> filterTumorCoverage(int count){
    if(count <= 10)
        return std::tuple<bool, std::string>(false, "low confidence call in plasma");
    if(count > 10 and count <= 50)
        return std::tuple<bool, std::string>(true, "low confidence VAF in plasma");
    else
        return std::tuple<bool, std::string>(true, "high confidence in plasma");
}

void importdbSNP(std::string dbSNPFilePath){
    std::ifstream dbSNPFile(dbSNPFilePath);

    std::string lineText;
    while ( std::getline(dbSNPFile, lineText)) {
        if(lineText[0] == '#')
            continue;

        std::vector<std::string> sp = split(lineText, "\t");

        if(sp[7].find("VC=SNV") == std::string::npos)
            continue;

        std::vector<std::string> spp = split(sp[2],",");

        std::string s = "ACGT";

        for(std::string& c : spp){
            if(c.size() > 1)
                continue;
            if(s.find(c) == std::string::npos){
                dbSNP["chr" + sp[0] + "-" + "A"].push_back(sp[1]);
                dbSNP["chr" + sp[0] + "-" + "C"].push_back(sp[1]);
                dbSNP["chr" + sp[0] + "-" + "G"].push_back(sp[1]);
                dbSNP["chr" + sp[0] + "-" + "T"].push_back(sp[1]);
            }
            else{
                dbSNP["chr" + sp[0] + "-" + c].push_back(sp[1]);
            }
        }
    }

    dbSNPFile.close();
}

bool filterSupportingCountByBinomialTest(double jenksEstimate, std::string basestring, char variantBase){
    int cnt = basestring.size();
    std::string basestringUpper = toUpper(basestring);
    int varCNT = std::count(basestringUpper.begin(), basestringUpper.end(), variantBase);
    //std::cout << varCNT << " " << cnt << " " << jenksEstimate/2 << std::endl; //?
    double p = binomialCDF(varCNT, cnt, jenksEstimate/2);
    if(p < P_VALUE_BINOMIAL_TEST_JENKS_ESTIMATE)
        return false;
    else
        return true;
}

std::tuple<double, int, int> finalEstimationWithJenks(std::vector<double> VAFList){
    if(VAFList.size() == 0) {
        //std::cout << "1st" << std::endl;
        return std::tuple<double, int, int>(0.0, 0, 0);
    }

    if(VAFList.size() < 20) {
        //std::cout << "2nd" << std::endl;
        return std::tuple<double, int, int>(sumVector(VAFList)/VAFList.size(), VAFList.size(), 0);
    }


    double tmp1 = 0.0;
    double tmp2 = 0.0;
    int maxIterLim = 10;
    int minIterLim = 3;
    int maxIterXLim = round(VAFList.size()/20.0);
    int maxLim = std::max(minIterLim, std::min(maxIterLim, maxIterXLim)) + 1;
    std::vector<std::vector<double>> groups;
    std::vector<double> breaks;

    //std::cout << minIterLim << " " << maxLim << std::endl;

    for( int ct = minIterLim; ct < maxLim; ct++){
        std::tuple<std::vector<std::vector<double>>, std::vector<double>> results = jenks(VAFList, ct);
        groups = std::get<0>(results);
        breaks = std::get<1>(results);
        //std::cout << groups.size() << std::endl;

        tmp2 = sumVector(groups[groups.size()-1])/groups[groups.size()-1].size();
        if(tmp1 == tmp2)
            break;
        else
            tmp1 = tmp2;
    }

    std::vector<double> include = {};
    int n = 0;
    for(int i = 0; i < groups.size(); i++){
        std::vector<double> group = groups[groups.size()-i-1];
        include.insert(include.end(), group.begin(), group.end());
        n += 1;
        if(include.size() >= VARIANT_COUNT_STOP_ADDING_IN_JENK_ESTIMATION)
            break;
    }

    if(include.size() > VARIANT_COUNT_FOR_JENK_ESTIMATION)
        include = {include.begin(), include.begin()+VARIANT_COUNT_FOR_JENK_ESTIMATION};

    //std::cout << "3rd" << include.size() << std::endl;

    return std::tuple<double, int, int>{sumVector(include)/include.size(), include.size(), n};
}
