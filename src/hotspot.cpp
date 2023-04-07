//
// Created by Colin Small on 7/6/21.
//

#include <Rcpp.h>
using namespace Rcpp;

#include "hotspot.h"

hotspot::hotspot(double priora, double priorb, double priorc, const std::string & chrom, const std::string & pos,
                 const std::string & basestring, const std::vector<double> & quallist,
                 const std::vector<double> & maplist, const std::string & basestringNormal,
                 const std::vector<double> & quallistNormal, const std::vector<double> & maplistNormal,
                 const std::string & basestringExtendedFrags, const std::vector<double> & quallistExtendedFrags,
                 const std::vector<double> & maplistExtendedFrags, const std::string & basestringNotCombined,
                 const std::vector<double> & quallistNotCombined, const std::vector<double> & maplistNotCombined,
                 char variantBase, int nread) : priora(priora), priorb(priorb), priorc(priorc), chrom(chrom), pos(pos), basestring(basestring), quallist(quallist),
                                     maplist(maplist), basestringNormal(basestringNormal),
                                     quallistNormal(quallistNormal), maplistNormal(maplistNormal),
                                     basestringExtendedFrags(basestringExtendedFrags),
                                     quallistExtendedFrags(quallistExtendedFrags),
                                     maplistExtendedFrags(maplistExtendedFrags),
                                     basestringNotCombined(basestringNotCombined),
                                     quallistNotCombined(quallistNotCombined), maplistNotCombined(maplistNotCombined),
                                     variantBase(variantBase), nread(nread) {

}
