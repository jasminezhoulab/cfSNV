//
// Created by Colin Small on 7/30/21.
//

#include <Rcpp.h>
using namespace Rcpp;

#include "filter_with_pileup.h"
#include "_helper.h"
#include "_filter.h"
#include "_probability.h"

#include <iostream>
#include <fstream>
#include <tuple>
#include <string>
#include <iomanip>

std::tuple<std::vector<std::string>, std::vector<bool>, std::string, bool, std::string, std::vector<bool>, std::string, bool, std::string, char, double> generateIntermediateResultOneLine(std::string line, double depth)
{
    char inVAF = 'F';
    std::vector<std::string> sp = split(line, "\t");

    for(int i = 0; i < sp.size(); i++)
        if(!sp[i].empty())
                sp[i] = sp[i].substr(1, sp[i].size()-2);

    std::string chrom = sp[0];
    std::string pos = sp[1];
    char variantBase = sp[9][0];
    double VAFMerge = std::stod(sp[14]);
    double VAF = std::stod(sp[11]);

    std::string basestring = sp[16];
    std::vector<double> quallist = stringToQual(sp[17]);
    std::vector<double> maplist = stringToQual(sp[18]);

    std::string basestringExtendedFrags = sp[22];
    std::vector<double> quallistExtendedFrags = stringToQual(sp[23]);
    std::vector<double> maplistExtendedFrags = stringToQual(sp[24]);

    std::string basestringNotCombined = sp[25];
    std::vector<double> quallistNotCombined = stringToQual(sp[26]);
    std::vector<double> maplistNotCombined = stringToQual(sp[27]);

    int cntNormal = sp[19].size();

    std::string basestringMerge = basestringExtendedFrags + basestringNotCombined;
    std::vector<double> quallistMerge;
    quallistMerge.reserve(quallistExtendedFrags.size() + quallistNotCombined.size());
    quallistMerge.insert(quallistMerge.end(), quallistExtendedFrags.begin(), quallistExtendedFrags.end());
    quallistMerge.insert(quallistMerge.end(), quallistNotCombined.begin(), quallistNotCombined.end());

    std::vector<double> maplistMerge;
    maplistMerge.reserve(maplistExtendedFrags.size() + maplistNotCombined.size());
    maplistMerge.insert(maplistMerge.end(), maplistExtendedFrags.begin(), maplistExtendedFrags.end());
    maplistMerge.insert(maplistMerge.end(), maplistNotCombined.begin(), maplistNotCombined.end());

    std::map<char, int> basecountExtendedFrags = countBase(basestringExtendedFrags);
    std::map<char, int> basecountNotCombined = countBase(basestringNotCombined);
    std::map<char, int> basecount = countBase(basestring);

    bool strandBiasMerge = filterStrandBiasMerge(basecountNotCombined, basecountExtendedFrags, variantBase);
    std::tuple<bool, bool> f = filterBothStrandAboveAverageMerge(basecountNotCombined, basecountExtendedFrags);
    bool bothStrandMerge = std::get<0>(f);
    bool aboveAverageMerge = std::get<1>(f);
    bool strandBiasUnmerge = filterStrandBiasUnmerge(basecount, variantBase);
    std::tuple<bool, bool> f2 = filterBothStrandAboveAverageUnmerge(basecount);
    bool bothStrandUnmerge = std::get<0>(f2);
    bool aboveAverageUnmerge = std::get<1>(f2);
    bool allMappingQualityUnmerge = filterAllMappingQuality(maplist);
    bool allMappingQualityMerge = filterAllMappingQuality(maplistMerge);
    bool variantMappingqualityUnmerge = filterVariantMappingQuality(maplist, basestring, variantBase);
    bool variantMappingqualityMerge = filterVariantMappingQuality(maplistMerge, basestringMerge, variantBase);
    bool variantBaseQualityUnmerge = filterVariantBaseQuality(quallist, basestring, variantBase);
    bool variantBaseQualityMerge = filterVariantBaseQuality(quallistMerge, basestringMerge, variantBase);
    bool supportingReadUnmerge = filterSupport(basestring, variantBase);
    bool supportingReadMerge = filterSupport(basestringMerge, variantBase);

    std::tuple<bool, std::string> VAFFilterUnmergeTuple = filterVAF(VAF, depth);
    bool VAFFilterUnmerge = std::get<0>(VAFFilterUnmergeTuple);
    std::string VAFDescriptionUnmerge = std::get<1>(VAFFilterUnmergeTuple);

    std::tuple<bool, std::string> VAFFilterMergeTuple = filterVAF(VAFMerge, depth);
    bool VAFFilterMerge = std::get<0>(VAFFilterMergeTuple);
    std::string VAFDescriptionMerge = std::get<1>(VAFFilterMergeTuple);

    int cnt = basestring.size();
    int cntMerge = basestringMerge.size();

    bool normalCoverage = filterNormalCoverage(cntNormal);

    std::tuple<bool, std::string> tumorCoverageUnmergeTuple = filterTumorCoverage(cnt);
    bool tumorCoverageUnmerge = std::get<0>(tumorCoverageUnmergeTuple);
    std::string tumorCoverageDescrptionUnmerge = std::get<1>(tumorCoverageUnmergeTuple);

    std::tuple<bool, std::string> tumorCoverageMergeTuple = filterTumorCoverage(cntMerge);
    bool tumorCoverageMerge = std::get<0>(tumorCoverageMergeTuple);
    std::string tumorCoverageDescrptionMerge = std::get<1>(tumorCoverageMergeTuple);

    if(normalCoverage and strandBiasMerge and bothStrandMerge and aboveAverageMerge and allMappingQualityMerge and variantMappingqualityMerge and variantBaseQualityMerge and supportingReadMerge and VAFFilterMerge
    and strandBiasUnmerge and bothStrandUnmerge and aboveAverageUnmerge and allMappingQualityUnmerge and variantMappingqualityUnmerge and variantBaseQualityUnmerge and supportingReadUnmerge and VAFFilterUnmerge
    and cnt > 50 and cntMerge > 50)
        inVAF = 'T';

    std::vector<std::string> outputStr;
    outputStr.reserve(28);
    outputStr.insert(outputStr.end(), sp.begin(), sp.begin() + 3);
    outputStr.push_back(sp[9]);
    outputStr.insert(outputStr.end(), sp.begin()+3, sp.begin()+7);
    outputStr.insert(outputStr.end(), sp.begin()+10, sp.begin()+13);
    outputStr.insert(outputStr.end(), sp.begin()+7, sp.begin()+9);
    outputStr.insert(outputStr.end(), sp.begin()+13, sp.begin()+16);
    outputStr.insert(outputStr.end(), sp.begin()+16, sp.begin()+28);

    std::vector<bool> outputMerge = {normalCoverage, strandBiasMerge, bothStrandMerge, aboveAverageMerge, allMappingQualityMerge, variantMappingqualityMerge, variantBaseQualityMerge, supportingReadMerge, VAFFilterMerge};
    std::vector<bool> outputUnmerge = {strandBiasUnmerge, bothStrandUnmerge, aboveAverageUnmerge, allMappingQualityUnmerge, variantMappingqualityUnmerge, variantBaseQualityUnmerge, supportingReadUnmerge, VAFFilterUnmerge};

    return std::tuple<std::vector<std::string>, std::vector<bool>, std::string, bool, std::string, std::vector<bool>, std::string, bool, std::string, char, double>{outputStr, outputMerge, VAFDescriptionMerge, tumorCoverageMerge, tumorCoverageDescrptionMerge, outputUnmerge, VAFDescriptionUnmerge, tumorCoverageUnmerge, tumorCoverageDescrptionUnmerge, inVAF, VAFMerge};
}

inline char boolToChar(bool b){
    return b ? 'T' : 'F';
}

std::vector<double> generateIntermediateResultWholeFile(std::string filename, std::string output, std::string database, double depth){
    importdbSNP(database);
    std::ifstream f(filename);
    std::ofstream g(output);

    std::vector<double> VAFList = {};

    std::string lineText;
    while (std::getline(f, lineText)) {
        std::tuple<std::vector<std::string>, std::vector<bool>, std::string, bool, std::string, std::vector<bool>, std::string, bool, std::string, char, double> lineTuple = generateIntermediateResultOneLine(lineText, depth);
        std::vector<std::string> outputStr = std::get<0>(lineTuple);
        std::vector<bool> outputMerge = std::get<1>(lineTuple);
        std::string VAFDescriptionMerge = std::get<2>(lineTuple);
        bool tumorCoverageMerge = std::get<3>(lineTuple);
        std::string tumorCoverageDescrptionMerge = std::get<4>(lineTuple);
        std::vector<bool> outputUnmerge = std::get<5>(lineTuple);
        std::string VAFDescriptionUnmerge = std::get<6>(lineTuple);
        bool tumorCoverageUnmerge = std::get<7>(lineTuple);
        std::string tumorCoverageDescrptionUnmerge = std::get<8>(lineTuple);
        char inVAF = std::get<9>(lineTuple);
        double VAFMerge = std::get<10>(lineTuple);

        if(inVAF == 'T')
            VAFList.push_back(VAFMerge);

        for(std::string& s : outputStr)
            g << s << '\t';

        for(bool b : outputMerge)
            g << boolToChar(b) << '\t';

        g << VAFDescriptionMerge << '\t';
        g << boolToChar(tumorCoverageMerge) << '\t';
        g << tumorCoverageDescrptionMerge << '\t';

        for(bool b : outputUnmerge)
            g << boolToChar(b) << '\t';

        g << VAFDescriptionUnmerge << '\t';
        g << boolToChar(tumorCoverageUnmerge) << '\t';
        g << tumorCoverageDescrptionUnmerge << '\t';
        g << '\n';
//        g << boolToChar(inVAF) << '\t';
//        g << VAFMerge << '\t' << '\n';
    }
    g.close();
    f.close();
    return VAFList;
}

const std::string headerIntermediateFile = "chromosome\tposition\treference\tvariant\t cnt_normal\t cnt_var_normal\t cnt_tumor\t cnt_var_tumor\tjoint_genotype_MAP\tVAF_tumor\tconfidence\t cnt_tumor_merge\t cnt_var_tumor_merge\tjoint_genotype_MAP_merge\tVAF_tumor_merge\tconfidence_merge\tbasestring\tqual_string\tmap_string\tbasestring_normal\tqual_string_normal\tmap_string_normal\tbasestring_ef\tqual_string_ef\tmap_string_ef\tbasestring_nc\tqual_string_nc\tmapstring_nc\tnormal_coverage\tstrand_bias_merge\tboth_strand_merge\tabove_average_merge\tall_mapping_quality_merge\tvariant_mapping_quality_merge\tvariant_base_quality_merge\tsupporting_frag_merge\tVAF_filter_merge\tVAF_descript_merge\ttumor_coverage_merge\tdescript_merge\tstrand_bias_unmerge\tboth_strand_unmerge\tabove_average_unmerge\tall_mapping_quality_unmerge\tvariant_mapping_quality_unmerge\tvariant_base_quality_unmerge\tsupporting_frag_unmerge\tVAF_filter_unmerge\tVAF_descript_unmerge\ttumor_coverage_unmerge\tdescript_unmerge";
const std::string headerRecordFile = "chromosome\tposition\treference\tvariant\t cnt_normal\t cnt_var_normal\t cnt_tumor\t cnt_var_tumor\tjoint_genotype_MAP\tVAF_tumor\tconfidence\tcnt_tumor_merge\tcnt_var_tumor_merge\tjoint_genotype_MAP_merge\tVAF_tumor_merge\tconfidence_merge\tbasestring\tqual_string\tmap_string\tbasestring_normal\tqual_string_normal\tmap_string_normal\tbasestring_ef\tqual_string_ef\tmap_string_ef\tbasestring_nc\tqual_string_nc\tmapstring_nc\tnormal_coverage\tstrand_bias_merge\tboth_strand_merge\tabove_average_merge\tall_mapping_quality_merge\tvariant_mapping_quality_merge\tvariant_base_quality_merge\tsupporting_frag_merge\tVAF_filter_merge\tVAF_descript_merge\ttumor_coverage_merge\tdescript_merge\tbinomial_test_merge\tstatus_merge\tstrand_bias_unmerge\tboth_strand_unmerge\tabove_average_unmerge\tall_mapping_quality_unmerge\tvariant_mapping_quality_unmerge\tvariant_base_quality_unmerge\tsupporting_frag_unmerge\tVAF_filter_unmerge\tVAF_descript_unmerge\ttumor_coverage_unmerge\tdescript_unmerge\tbinomial_test_unmerge\tstatus_unmerge\tdecision";
const std::string headerOutputFile = "chromosome\tposition\treference\tvariant\tcnt_normal\tcnt_var_normal\tcnt_tumor\tcnt_var_tumor\tjoint_genotype_MAP\tVAF_tumor\tconfidence\tcnt_tumor_merge\t cnt_var_tumor_merge\tjoint_genotype_MAP_merge\tVAF_tumor_merge\tconfidence_merge";


std::vector<std::string> generateRecordOneLine(std::string line, double jenks_estimate){
    std::vector<std::string> sp = split(line, "\t");
    sp.pop_back();
    std::string chrom = sp[0];
    std::string pos = sp[1];
    std::string ref = sp[2];
    char variantBase = sp[3][0];
    std::string basestring = sp[16];
    std::string basestringExtendedFrags = sp[22];
    std::string basestringNotCombined = sp[25];
    std::string basestringMerge = basestringExtendedFrags + basestringNotCombined;
    bool binomTestUnmerge = filterSupportingCountByBinomialTest(jenks_estimate, basestring, variantBase);
    bool binomTestMerge = filterSupportingCountByBinomialTest(jenks_estimate, basestringMerge, variantBase);
    sp.insert(sp.begin()+40, std::string(1, boolToChar(binomTestMerge)));
    sp.emplace_back(std::string(1, boolToChar(binomTestUnmerge)));
    return sp;
}

std::tuple<std::vector<std::string>, std::string, std::string, std::string> decideOutputStatusOneLine(std::vector<std::string> record){
    std::vector<std::string> merge = {};
    merge.insert(merge.end(), record.begin()+29, record.begin()+37);
    merge.emplace_back(record[38]);
    merge.emplace_back(record[40]);

    std::vector<std::string> unmerge = {};
    unmerge.insert(unmerge.end(), record.begin()+41, record.begin()+49);
    unmerge.emplace_back(record[50]);
    unmerge.emplace_back(record[52]);

    std::vector<std::string> output = {record.begin(), record.begin()+16};

    std::string mergeStatus;
    std::string unmergeStatus;

    if(record[28] == "T" and record[41] == "T" and merge[0] == "T" and merge[2] == "T" and merge[3] == "T" and merge[4] == "T" and merge[5] == "T" and merge[6] == "T" and merge[7] == "T" and merge[8] == "T" and merge[9] == "T")
        mergeStatus = "pass";
    else if(record[28] == "T" and record[41] == "T" and merge[0] == "T" and merge[2] == "T" and merge[3] == "T" and merge[4] == "T" and merge[5] == "T" and merge[9] == "T")
        mergeStatus = "hold";
    else
        mergeStatus = "reject";

    if(record[28] == "T" and record[41] == "T" and unmerge[0] == "T" and unmerge[2] == "T" and unmerge[3] == "T" and unmerge[4] == "T" and unmerge[5] == "T" and unmerge[6] == "T" and unmerge[7] == "T" and unmerge[8] == "T" and unmerge[9] == "T")
        unmergeStatus = "pass";
    else if(record[28] == "T" and record[41] == "T" and unmerge[0] == "T" and unmerge[2] == "T" and unmerge[3] == "T" and unmerge[4] == "T" and unmerge[5] == "T" and unmerge[9] == "T")
        unmergeStatus = "hold";
    else
        unmergeStatus = "reject";

    if(mergeStatus == "pass" and unmergeStatus == "pass")
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "pass");
    else if( mergeStatus == "pass" and unmergeStatus == "hold")
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "hold");
    else if( mergeStatus == "hold" and unmergeStatus == "pass")
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "hold");
    else if( mergeStatus == "hold" and unmergeStatus == "hold")
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "hold");
    else if( mergeStatus == "pass" and unmergeStatus == "reject")
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "hold");
    else if( mergeStatus == "reject" and unmergeStatus == "pass")
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "hold");
    else if( mergeStatus == "reject" and unmergeStatus == "reject")
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "reject");
    else if( mergeStatus == "hold" and unmergeStatus == "reject")
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "hold");
    else
        return std::tuple<std::vector<std::string>, std::string, std::string, std::string>( output, mergeStatus, unmergeStatus, "hold");
}

int generateRecordAndOutputWholeFile(std::string fintermediate, std::string foutputPass, std::string foutputCheck, std::string frecord, double jenksEstimate){
    std::ifstream f(fintermediate);
    std::ofstream r(frecord);
    std::ofstream g(foutputPass);
    std::ofstream h(foutputCheck);

    r << headerRecordFile << '\n';
    int n = 0;

    std::string lineText;
    while ( std::getline(f, lineText)) {
        std::vector<std::string> record = generateRecordOneLine(lineText, jenksEstimate);
        std::tuple<std::vector<std::string>, std::string, std::string, std::string> outputStatusTuple = decideOutputStatusOneLine(record);
        std::vector<std::string> output = std::get<0>(outputStatusTuple);
        std::string mergeStatus = std::get<1>(outputStatusTuple);
        std::string unmergeStatus = std::get<2>(outputStatusTuple);
        std::string status = std::get<3>(outputStatusTuple);

        if(status == "pass") {
            n += 1;
            for(std::string& s : output){
                g << s << '\t';
            }
            g << '\n';
        }

        if(status == "hold") {
            n += 1;
            for(std::string& s : output){
                h << s << '\t';
            }
            h << '\n';
        }

        record.insert(record.begin()+41, mergeStatus);
        record.push_back(unmergeStatus);
        record.push_back(status);

        for(std::string& s : record){
            r << '"' << s << '"' << '\t';
        }
        r << '\n';
    }

    f.close();
    r.close();
    g.close();
    h.close();

    return n;
}

// [[Rcpp::export]]
int filter_with_pileup_main(std::string fInput, std::string fIntermediate, std::string fOutputPass,
                            std::string fOutputCheck, std::string fRecord, std::string fTF,
                            double mergedVAFThreshold, std::string SNPDatabase, double depth){
    // std::string fInput = argv[1];
    // std::string fIntermediate = argv[2];
    // std::string fOutputPass  = argv[3];
    // std::string fOutputCheck = argv[4];
    // std::string fRecord = argv[5];
    // std::string fTF = argv[6];
    // std::string mergedVAFThreshold = argv[7];
    // std::string SNPDatabase = argv[8];
    // double depth = std::stod(argv[9]);

    std::vector<double> VAFList = generateIntermediateResultWholeFile(fInput, fIntermediate, SNPDatabase, depth);
    auto t = finalEstimationWithJenks(VAFList);
    double jenksEstimate = std::get<0>(t);
    int includeNumber = std::get<1>(t);
    int includeGroup = std::get<2>(t);
    int n = generateRecordAndOutputWholeFile(fIntermediate, fOutputPass, fOutputCheck, fRecord, jenksEstimate);

    std::ofstream TF(fTF);
    TF << n << '\t' << std::setprecision(13) << std::fixed << jenksEstimate << '\t' << std::setprecision(4) << std::fixed << jenksEstimate/2.0 << '\t' << includeNumber << '\t' << includeGroup << '\t' << mergedVAFThreshold << '\n';
    TF.close();

    return 0;
}
