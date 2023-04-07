import math
from decimal import *
from numpy import median
from numpy import std
import numpy as np
import re
from scipy import stats

# GLOBAL VARIABLES

BASEQUAL_THRESHOLD = 0.05 # _probability
MAPQUAL_THRESHOLD = 0.3 # _probability
COUNT_VAR_TUMOR_HIGHQUAL = 3 # py1


#####################################################


DEPTH_FOR_GERMLINE = 30	# py2
GERMLINE_VAF_THRESHOLD_IN_NORMAL = 0.01 # py2
GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER = 2 # py2
GERMLINE_VARIANT_COUNT = 2 ### BAMSurgeon pooled 2 ### WES 1 # py2
MAX_NUMBER_OF_HOTSPOT = 50 # py2
AVG_MAPQUAL_FOR_ESTIMATION = 0.1 # py2
ZERO_MAPQUAL_COUNT_FOR_ESTIMATION = 5 #py2

GRIDWIDTH = 0.001 # py2
MAXSEARCH = 400 # py2
MAXSEARCH_WITH_NORMAL = 1000 # py2
DEPTH_FOR_ESTIMATION = 80 # py2


STRAND_BIAS_RATIO_VARIANT_TO_ALL = 5.0 # _filter
STRAND_BIAS_BINOMIAL_PROB = 0.05 # _filter
THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF = 0.8 # _filter
COUNT_FOR_STRONG_EVIDENCE = 3 # _filter

TRIALLELE_VAF_RATIO = 0.5 ### BAMSurgeon pooled cfDNA 1 ### BAMSurgeon pooled 0.5 ### WES 0.5 #_filter
TRIALLELE_VAF = 0.02 #_filter
TRIALLELE_COUNT = 4 ### BAMSurgeon pooled 4 ### WES 3 #_filter

## fixed
HOTSPOT = []
HOTSPOT_PRIOR = []

## prior probability for joint genotypes
## order: AA/AA, AA/AB, AA/BB,
##                AB/AA, AB/AB, AB/BB,
##                BB/AA, BB/AB, BB/BB.
EST_PRIOR = dict([])
EST_PRIOR["AA/AA"] = 10.0**6.0
EST_PRIOR["AA/AB"] = 10.0**4.0
EST_PRIOR["AA/BB"] = 10.0**2.0
EST_PRIOR["AB/AA"] = 10.0**2.0
EST_PRIOR["AB/AB"] = 10.0**3.0
EST_PRIOR["AB/BB"] = 0.0001
EST_PRIOR["BB/AA"] = 0.0001
EST_PRIOR["BB/AB"] = 0.0001
EST_PRIOR["BB/BB"] = 0.0001



##############################################################################


NORMAL_COUNT_ALT = 7 # py3
#NORMAL_COUNT_VAR = 2 # py3
SOMATIC_VAF_THRESHOLD_IN_NORMAL = 0.12 # py3
NORMAL_COUNT_BINOM = 0.2 ### BAMSurgeon pooled 0.2 ### WES 0.05 # py3
NORMAL_COUNT_FOR_VAF_THRESHOLD = 50 # py3
DEPTH_FOR_DETECTION_NORMAL = 4 # py3
DEPTH_FOR_DETECTION_TUMOR = 4 # py3

PRIOR = dict([])
PRIOR["AA/AA"] = 10.0**6.0
PRIOR["AA/AB"] = 10.0**4.0
PRIOR["AA/BB"] = 10.0**2.0
PRIOR["AB/AA"] = 10.0**2.0
PRIOR["AB/AB"] = 10.0**4.0
PRIOR["AB/BB"] = 10.0**2.0
PRIOR["BB/AA"] = 1.0
PRIOR["BB/AB"] = 1.0
PRIOR["BB/BB"] = 10.0**4.0





#############################################################################



VAF_FOR_STRONG_EVIDENCE = 0.0005 ### BAMSurgeon pooled 0.0005 ### WES 0.015
LOW_MAPQUAL_COUNT_FOR_DETECTION_VARIANT = 3
DEPTH_ACTIVATE_LOW_MAPQUAL_COUNT = 6
LOW_MAPQUAL_COUNT_FOR_DETECTION_ALL = 10
DEPTH_ACTIVATE_LOW_MAPQUAL_COUNT_ALL = 20
LOW_MAPQUAL = 0.3
OK_MAPQUAL = 0.1
OK_BASEQUAL = 0.005
THRESHOLD_OK_BASEQUAL_COUNT = 3
THRESHOLD_OK_MAPQUAL_COUNT = 3
THRESHOLD_OK_MAPQUAL_COUNT_ALL = 20
AVG_MAPQUAL_FOR_DETECTION = 0.1
LOW_MAPQUAL_FRAC_FOR_DETECTION = 0.4
THRESHOLD_VARIANT_QUALITY_TO_REFERENCE_QUALITY = 7.0
VARIANT_COUNT_FOR_JENK_ESTIMATION = 30 ### Zo 50 ### Viktor 30
P_VALUE_BINOMIAL_TEST_JENKS_ESTIMATE = 0.1
VARIANT_COUNT_STOP_ADDING_IN_JENK_ESTIMATION = 20 ### Viktor 20 ### Zo 50
# post_filter2: find_all_occurrence_in_string
CIGAR_CHARACTER = "MSIDHNPX="
# post_filter2: if_indel_near
DISTANCE_FROM_INDEL = 4
# post_filter2: filter_quality_and_indel
NUMBER_INDEL_NEARBY = 3
THRESHOLD_ADJACENT_ALL_T_STATISTICS = 10
THRESHOLD_ADJACENT_VARIANT_T_STATISTICS = 10



#############################################################################



