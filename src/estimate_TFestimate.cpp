//
// Created by Colin Small on 7/6/21.
//

#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <iomanip>

#include "estimate_TFestimate.h"
#include "_helper.h"
#include "_parameter.h"
#include "_probability.h"
#include "_filter.h"

//double MERGED_VAF_THRESHOLD;

// ZERO_MAPQUAL_COUNT_FOR_ESTIMATION;
const int MAX_NUMBER_OF_HOTPOST = 50;

void importPredefinedHotspotPrior(std::fstream file) {

    std::string line;
    std::vector<std::string> splitVector;
    while ( getline(file, line)) {
        splitVector = split(line, "\t");
        HOTSPOT_PRIOR.emplace_back(splitVector[0], std::stoi(splitVector[1]));
    }
}

void selectHotSpot(std::string chrom, std::string pos, std::string basestring, std::vector<double> quallist,
                   std::vector<double> maplist, std::string basestringNormal, std::vector<double> quallistNormal,
                   std::vector<double> maplistNormal, std::string basestringExtendedFrags,
                   std::vector<double> quallistExtendedFrags, std::vector<double> maplistExtendedFrags,
                   std::string basestringNotCombined, std::vector<double> quallistNotCombined,
                   std::vector<double> maplistNotCombined, std::map<char, int> basecountNotCombined,
                   std::map<char, int> basecountExtendedFrags, char variantBase, int nread, double MERGED_VAF_THRESHOLD) {

    if ( nread < DEPTH_FOR_ESTIMATION )
        return;

    std::string basestringNormalUpper = toUpper(basestringNormal);

    double basestringNormalVariantCount = std::count(basestringNormalUpper.begin(), basestringNormalUpper.end(),
                                                     variantBase);

    if ( basestringNormalVariantCount > GERMLINE_VARIANT_COUNT )
        return;

    std::map<char, int> basecount = countBase(basestring);
    std::string basestringMerge = basestringExtendedFrags + basestringNotCombined;

    std::vector<double> quallistMerge;
    quallistMerge.reserve(quallistExtendedFrags.size() + quallistNotCombined.size());
    quallistMerge.insert(quallistMerge.end(), quallistExtendedFrags.begin(), quallistExtendedFrags.end());
    quallistMerge.insert(quallistMerge.end(), quallistNotCombined.begin(), quallistNotCombined.end());

    std::vector<double> maplistMerge;
    maplistMerge.reserve(maplistExtendedFrags.size() + maplistNotCombined.size());
    maplistMerge.insert(maplistMerge.end(), maplistExtendedFrags.begin(), maplistExtendedFrags.end());
    maplistMerge.insert(maplistMerge.end(), maplistNotCombined.begin(), maplistNotCombined.end());

    // Require DEPTH_FOR_ESTIMATION fragments
    if ( basestringMerge.size() < DEPTH_FOR_ESTIMATION )
        return;
    // Exclude loci without enough germline information
    if ( basestringNormal.size() < DEPTH_FOR_GERMLINE )
        return;
    // Require high-quality variant alleles
    if ( basestring.size() == 0 || basestringMerge.size() == 0 || basestringNormal.size() == 0 )
        return;

    std::string basestringUpper = toUpper(basestring);
    int countNonRefRead = std::count(basestringUpper.begin(), basestringUpper.end(), variantBase);
    double meanQual = sumVector(quallist) / (double) quallist.size();

    if ((double) countNonRefRead < meanQual * (double) basestring.size())
        return;

    std::string baseStringMergeUpper = toUpper(basestringMerge);
    int countNonRefReadMerge = std::count(baseStringMergeUpper.begin(), baseStringMergeUpper.end(), variantBase);

    double meanQuallistMerge = sumVector(quallistMerge) / (double) quallistMerge.size();

    if ((double) countNonRefReadMerge < meanQuallistMerge * (double) basestringMerge.size())
        return;

    // Merged VAF lower than threshold
    // TODO: This corresponds to line 56 in p2.estimate.TFestimate.py, make sure that cnt_var_read_merge in that file
    // TODO: is the same as cnt_non_ref_read_merge and not a typo !!! nonono
    if ((double) countNonRefReadMerge / (double) basestringMerge.size() > MERGED_VAF_THRESHOLD )
        return;

    // Exclude possible germline variants
    if ((double) basestringNormalVariantCount / (double) basestringNormal.size() > GERMLINE_VAF_THRESHOLD_IN_NORMAL )
        return;

    // Require high maping quality in both tumor and normal
    int maplistZeroCount = std::count(maplist.begin(), maplist.end(), 0);
    int maplistNormalZeroCount = std::count(maplistNormal.begin(), maplistNormal.end(), 0);
    int maplistMergeZeroCount = std::count(maplistMerge.begin(), maplistMerge.end(), 0);
    if ( maplistZeroCount > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION ||
         maplistNormalZeroCount > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION ||
         maplistMergeZeroCount > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION )
        return;

    std::vector<double> maplistGreaterThanZero = accessMultipleIndicesVector(maplist,
                                                                             getIndicesWhereGreaterThan(maplist, 0.0));
    std::vector<double> maplistNormalGreaterThanZero = accessMultipleIndicesVector(maplistNormal,
                                                                                   getIndicesWhereGreaterThan(
                                                                                           maplistNormal, 0.0));
    std::vector<double> maplistMergeGreaterThanZero = accessMultipleIndicesVector(maplistMerge,
                                                                                  getIndicesWhereGreaterThan(
                                                                                          maplistMerge, 0.0));
    double maplistMean = meanVector(maplistGreaterThanZero);
    double maplistNormalMean = meanVector(maplistNormalGreaterThanZero);
    double maplistMergeMean = meanVector(maplistMergeGreaterThanZero);
    if ( maplistMean > AVG_MAPQUAL_FOR_ESTIMATION ||
         maplistNormalMean > AVG_MAPQUAL_FOR_ESTIMATION ||
         maplistMergeMean > AVG_MAPQUAL_FOR_ESTIMATION )
        return;

    // Exclude strand bias
    bool strandBiasMerge = filterStrandBiasMerge(basecountNotCombined, basecountExtendedFrags, variantBase);

    // Require enough percentage in non-reference reads
    std::tuple<bool, bool> bothObservedAboveAverageMerge = filterBothStrandAboveAverageMerge(basecountNotCombined,
                                                                                             basecountExtendedFrags);
    bool bothObservedMerge = std::get<0>(bothObservedAboveAverageMerge);
    bool aboveAverageMerge = std::get<1>(bothObservedAboveAverageMerge);

    bool supportFragCountMerge = filterSupport(basestringMerge, variantBase);
    bool strandBiasUnmerge = filterStrandBiasUnmerge(basecount, variantBase);

    // Require enough percentage in non-reference reads
    std::tuple<bool, bool> bothObservedAboveAverageUnmerge = filterBothStrandAboveAverageUnmerge(basecount);
    bool bothObservedUnmerge = std::get<0>(bothObservedAboveAverageUnmerge);
    bool aboveAverageUnmerge = std::get<1>(bothObservedAboveAverageUnmerge);

    bool supportReadCountUnmerge = filterSupport(basestring, variantBase);

    if ( !strandBiasMerge or !bothObservedMerge or !aboveAverageMerge or !supportFragCountMerge )
        return;
    if ( !strandBiasUnmerge or !bothObservedUnmerge or !aboveAverageUnmerge or !supportReadCountUnmerge)
        return;

    HOTSPOT.emplace_back(hotspot(meanQuallistMerge, basestringNormalVariantCount,
                         (double) countNonRefReadMerge / (double) basestringMerge.size(),
                         chrom, pos, basestring, quallist, maplist, basestringNormal, quallistNormal, maplistNormal,
                         basestringExtendedFrags, quallistExtendedFrags, maplistExtendedFrags, basestringNotCombined,
                         quallistNotCombined, maplistNotCombined, variantBase, nread));
}

// Unused in original python script
//void selectHotspotFromPredefinedHotspotPrior(std::string_view chrom, int pos, std::string_view basestring,
//                                             std::vector<double> & quallist, std::vector<double> & maplist,
//                                             std::string_view basestringNormal, std::vector<double> & quallistNormal,
//                                             std::vector<double> & maplistNormal,
//                                             std::string_view basestringExtendedFrags,
//                                             std::vector<double> & quallistExtendedFrags,
//                                             std::vector<double> & maplistExtendedFrags,
//                                             std::string_view basestringNotCombined,
//                                             std::vector<double> & quallistNotCombined,
//                                             std::vector<double> & maplistNotCombined,
//                                             std::map<char, int> & basecountNotCombined,
//                                             std::map<char, int> & basecountExtendedFrags, char variantBase) {
//    if( vecContains(HOTSPOT_PRIOR, std::tuple<std::string, int>(chrom, pos))){
//        std::cout << "(" << chrom << ", " << pos << ")";
//        selectHotSpot(chrom, pos, basestring, quallist, maplist, basestringNormal, quallistNormal, maplistNormal,  basestringExtendedFrags, quallistExtendedFrags, maplistExtendedFrags, basestringNotCombined, quallistNotCombined, maplistNotCombined, basecountNotCombined, basecountExtendedFrags, variantBase);
//    }
//}

double calculateTumorFractionLikelihood(double tumor_fraction, std::string basestring_merge,
                                        std::vector<double> & quallist_merge, std::vector<double> & maplist_merge,
                                        std::string basestring_normal, std::vector<double> & quallist_normal,
                                        std::vector<double> & maplist_normal, char variant_base) {

    double likelihood = 0.0;
    std::map<std::string, double>::iterator it;
    std::string index;
    double loglikelihood;
    for ( it = PRIOR.begin() ; it != PRIOR.end() ; it++ ) {
        index = it->first;
        double l1 = calculateJointGenotypeTumorFractionLoglikelihood(tumor_fraction, basestring_merge, quallist_merge,
                                                              maplist_merge, index, variant_base);
        double l2 = calculateJointGenotypeTumorFractionLoglikelihood(0.0, basestring_normal, quallist_normal,
                                                                     maplist_normal, index, variant_base);
        loglikelihood = l1 + l2;
        likelihood += exp(loglikelihood + std::log(EST_PRIOR.find(index)->second));
    }
    return likelihood;
}

std::tuple<double, std::vector<double>> estimateTumorFraction(std::vector<hotspot> HOTSPOT) {

    std::vector<double> ratio(MAXSEARCH_WITH_NORMAL, 0);
    for ( int i = 0 ; i < MAXSEARCH_WITH_NORMAL ; i++ ) {
        ratio[i] = (double)GRIDWIDTH * i;
    }

    std::vector<double> tumorFractionLikelihood(MAXSEARCH_WITH_NORMAL, 0);
    for ( int ratioind = 0 ; ratioind < MAXSEARCH_WITH_NORMAL ; ratioind++ ) {
        for ( hotspot spot : HOTSPOT ) {
            std::string basestringMerge = spot.basestringExtendedFrags + spot.basestringNotCombined;
            std::vector<double> quallistMerge = spot.quallistNotCombined;
            quallistMerge.insert(quallistMerge.begin(), spot.quallistExtendedFrags.begin(),
                                 spot.quallistExtendedFrags.end());
            std::vector<double> maplistMerge = spot.maplistNotCombined;
            maplistMerge.insert(maplistMerge.begin(), spot.maplistExtendedFrags.begin(), spot.maplistExtendedFrags.end());

            double likelihood = calculateTumorFractionLikelihood(ratio[ratioind], basestringMerge, quallistMerge,
                                                                 maplistMerge, spot.basestringNormal,
                                                                 spot.quallistNormal, spot.maplistNormal,
                                                                 spot.variantBase);
            tumorFractionLikelihood[ratioind] += std::log(likelihood) / (spot.basestringExtendedFrags.size() +
                                                                         spot.basestringNotCombined.size() +
                                                                         spot.basestringNormal.size());
        }
    }

    int estind = std::max_element(tumorFractionLikelihood.begin(), tumorFractionLikelihood.end()) -
                 tumorFractionLikelihood.begin();
    double est = ratio[estind];
    return std::tuple<double, std::vector<double>>(est, tumorFractionLikelihood);
}

bool hotsportSortHelper(hotspot a, hotspot b){

    if(a.priora == b.priora){
        if(a.priorc == b.priorc){
            return a.priorb < b.priorb;
        }
        else
            return a.priorc > b.priorc;
    }
    else
        return a.priora < b.priora;
}

bool hotsportSort0(hotspot a, hotspot b) {
    return a.priora < b.priora;
}

bool hotsportSort1(hotspot a, hotspot b) {
    return  a.priorb <  b.priorb;
}

bool hotsportSort2(hotspot a, hotspot b) {
    return  b.priorc <  a.priorc;
}


// [[Rcpp::export]]
int estimate_TFestimate_main(std::string filename, double MERGED_VAF_THRESHOLD, std::string file_prefix,
                             double depth, std::string VAF_output, std::string estimate_output) {
    // std::cout << filename << std::endl;
    // std::cout << MERGED_VAF_THRESHOLD << std::endl;
    // std::cout << file_prefix << std::endl;
    // std::cout << depth << std::endl;
    // std::cout << VAF_output << std::endl;
    // std::cout << estimate_output << std::endl;

    HOTSPOT = {};
    HOTSPOT_PRIOR = {};

    GERMLINE_VARIANT_COUNT = getGermlineVariantCountThreshold(depth);

    std::ifstream file(filename);

    std::string lineText;
    std::string pos;

    while ( std::getline(file, lineText)) {
        std::vector<std::string> splitString = split(lineText, "\t");
        std::string chrom = splitString[0];

        if ( chrom.find('X') != std::string::npos or
             chrom.find('M') != std::string::npos or
             chrom.find('Y') != std::string::npos )
            continue;

        pos = splitString[1];
        int nread = splitString[19].size();
        std::string basestringAll = splitString[4];
        std::string basestring = splitString[19];

        if ( basestring.size() == 0 )
            continue;

        std::vector<double> quallist = stringToQual(splitString[20]);
        std::vector<double> maplist = stringToQual(splitString[21]);

        if ( splitString.size() < 31 )
            continue;

        std::string basestring_normal_all = splitString[8];
        std::string basestring_normal = splitString[22];
        std::vector<double> quallist_normal = stringToQual(splitString[23]);
        std::vector<double> maplist_normal = stringToQual(splitString[24]);
        std::string basestring_extendedFrags = splitString[25];
        std::vector<double> quallist_extendedFrags = stringToQual(splitString[26]);
        std::vector<double> maplist_extendedFrags = stringToQual(splitString[27]);
        std::string basestring_notCombined = splitString[28];
        std::vector<double> quallist_notCombined = stringToQual(splitString[29]);
        std::vector<double> maplist_notCombined = stringToQual(splitString[30]);
        std::map<char, int> basecount_notCombined_all = countBase(splitString[12]);
        std::map<char, int> basecount_extendedFrags_all = countBase(splitString[16]);
        char variantBaseAll = findMajorVariant(basecount_notCombined_all, basecount_extendedFrags_all);
        std::map<char, int> basecountNotCombined = countBase(basestring_notCombined);
        std::map<char, int> basecountExtendedFrags = countBase(basestring_extendedFrags);
        char variantBase = findMajorVariant(basecountNotCombined, basecountExtendedFrags);
        std::map<char, int> basecountTumorAll = countBase(basestringAll);
        std::map<char, int> basecountTumor = countBase(basestring);
        bool trialleleTumorAll = filterTriallelicPosition(basecountTumorAll, variantBaseAll, depth);
        bool trialleleTumor = filterTriallelicPosition(basecountTumor, variantBaseAll, depth);

        if ( variantBase != variantBaseAll )
            continue;
        if ( std::count(basestring.begin(), basestring.end(), std::toupper(variantBase)) +
             std::count(basestring.begin(), basestring.end(), std::tolower(variantBase)) == 0 )
            continue;
        if ( !trialleleTumorAll or !trialleleTumor )
            continue;

        std::string basestringNormallAllUpper = toUpper(basestring_normal_all);
        if ( std::count(basestringNormallAllUpper.begin(), basestringNormallAllUpper.end(), variantBase) >
             GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER )
            continue;
        if ( std::count(basestringNormallAllUpper.begin(), basestringNormallAllUpper.end(), variantBaseAll) >
             GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER )
            continue;
        selectHotSpot(chrom, pos, basestring, quallist, maplist, basestring_normal, quallist_normal, maplist_normal,
                      basestring_extendedFrags, quallist_extendedFrags, maplist_extendedFrags, basestring_notCombined,
                      quallist_notCombined, maplist_notCombined, basecountNotCombined, basecountExtendedFrags,
                      variantBase, nread, MERGED_VAF_THRESHOLD);
    }
    file.close();

    if ( HOTSPOT.size() > MAX_NUMBER_OF_HOTPOST ) {
        std::stable_sort(HOTSPOT.begin(), HOTSPOT.end(), hotsportSort0);
        std::stable_sort(HOTSPOT.begin(), HOTSPOT.end(), hotsportSort2);
        std::stable_sort(HOTSPOT.begin(), HOTSPOT.end(), hotsportSort1);
        HOTSPOT = {HOTSPOT.begin(), HOTSPOT.begin() + MAX_NUMBER_OF_HOTPOST};
    }

    std::tuple<double, std::vector<double>> estTumorFractionLikelihood = estimateTumorFraction(HOTSPOT);
    double est = std::get<0>(estTumorFractionLikelihood);
    std::vector<double> tumorFractionLikelihood = std::get<1>(estTumorFractionLikelihood);

    std::vector<double> VAF;
    VAF.reserve(HOTSPOT.size());

    // std::vector<hotspot> h = HOTSPOT;

    for ( hotspot & hs : HOTSPOT ) {
        std::string basestringCombined = hs.basestringNotCombined + hs.basestringExtendedFrags;
        int basecount = std::count(basestringCombined.begin(), basestringCombined.end(), std::toupper(hs.variantBase)) +
                        std::count(basestringCombined.begin(), basestringCombined.end(), std::tolower(hs.variantBase));
        double vaf = (double) basecount / basestringCombined.size();
        VAF.push_back(vaf);
    }
    std::sort(VAF.begin(), VAF.end());


    std::ofstream vafOutputFile(VAF_output);
    for ( double & vaf : VAF )
        vafOutputFile << std::setprecision(12) << std::fixed << vaf << '\n';
    vafOutputFile.close();

    std::ofstream estimateOutputFile(estimate_output);
    estimateOutputFile << est << '\n';
    estimateOutputFile << "====================" << '\n';
    for ( double & tfl : tumorFractionLikelihood )
        estimateOutputFile << std::setprecision(10) << std::fixed << tfl << '\n';
    estimateOutputFile.close();
    return 0;
}
