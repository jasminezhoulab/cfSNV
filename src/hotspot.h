//
// Created by Colin Small on 7/6/21.
//

#ifndef CFSNV_HOTSPOT_H
#define CFSNV_HOTSPOT_H

#include <Rcpp.h>
using namespace Rcpp;

#include <tuple>
#include <string>
#include <vector>

class hotspot{
public:
    hotspot(double priora, double priorb, double priorc, const std::string & chrom, const std::string & pos,
            const std::string & basestring, const std::vector<double> & quallist, const std::vector<double> & maplist,
            const std::string & basestringNormal, const std::vector<double> & quallistNormal,
            const std::vector<double> & maplistNormal, const std::string & basestringExtendedFrags,
            const std::vector<double> & quallistExtendedFrags, const std::vector<double> & maplistExtendedFrags,
            const std::string & basestringNotCombined, const std::vector<double> & quallistNotCombined,
            const std::vector<double> & maplistNotCombined, char variantBase, int nread);

    double priora; //0
    double priorb; //1
    double priorc; //2
    std::string chrom; //3
    std::string pos; //4
    std::string basestring; //5
    std::vector<double> quallist; //6
    std::vector<double> maplist; //7
    std::string basestringNormal; //8
    std::vector<double> quallistNormal; //9
    std::vector<double> maplistNormal; //10
    std::string basestringExtendedFrags; //11
    std::vector<double> quallistExtendedFrags; //12
    std::vector<double> maplistExtendedFrags; //13
    std::string basestringNotCombined; //14
    std::vector<double> quallistNotCombined; //15
    std::vector<double> maplistNotCombined; //16
    char variantBase; //17
    int nread;
};


#endif //CFSNV_HOTSPOT_H
