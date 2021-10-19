from parameter import *
from _jenks import *
from _probability import *
import numpy as np
from collections import defaultdict
from scipy.stats import binom


def find_major_variant(basecount_nc, basecount_ef):
# input: converted basestring
# output: variant with the max count of observations
# 	haven't thought of the case where the max count corresponds to a variant with strong strand bias
	alternative_nucleotide = ['A', 'C', 'T', 'G']
	count = []
	for nucleotide in alternative_nucleotide:
		count.append( basecount_nc[nucleotide] + basecount_nc[nucleotide.lower()] + 2*basecount_ef[nucleotide] + 2*basecount_ef[nucleotide.lower()] )
	max_count = max(count)
	variant_base = alternative_nucleotide[count.index(max_count)]
	return variant_base


def get_germline_variant_count_threshold(depth):
	return int(math.log10(depth))-1



def get_normal_count_binom_threshold(depth):
	if depth < 1000:
		return 0.05
	else:
		return 0.2


def filter_triallelic_position(basecount, var, depth):
	var_count = basecount[var] + basecount[var.lower()]
	total_count = sum(basecount.values())
	if total_count == 0:
		print(total_count, var_count, basecount, var)
	var_VAF = float(var_count)/float(total_count)
	TRIALLELE_COUNT = int(math.log10(depth)) + 2
	other_count = 0
	other_VAF = 0
	for i in "ACGT":
		if i == var:
			continue
		tmp_count = basecount[i] + basecount[i.lower()]
		if tmp_count > other_count:
			other_count = tmp_count
	other_VAF = float(other_count)/float(total_count)
	if var_count == 0 or other_VAF/var_VAF > TRIALLELE_VAF_RATIO:
#		print basecount, var
		return "F"
	if other_VAF > TRIALLELE_VAF and total_count >= 80:
		return "F"
	if other_count > TRIALLELE_COUNT and total_count <= 140:
		return "F"
	return "T"

def filter_strand_bias_merge(basecount_nc, basecount_ef, var):
# input: base counts
# output: if the locus passes the filter
		# the percentage of reads from forward strand
# related global variables:
	# TOTAL_STRAND_BIAS_RATIO_HIGH
	# TOTAL_STRAND_BIAS_RATIO_LOW
	# print basecount_nc, basecount_ef
	forward = basecount_nc['R'] + basecount_nc['A'] + basecount_nc['C'] + basecount_nc['T'] + basecount_nc['G'] + (basecount_ef['r'] + basecount_ef['a'] + basecount_ef['c'] + basecount_ef['t'] + basecount_ef['g'] + basecount_ef['R'] + basecount_ef['A'] + basecount_ef['C'] + basecount_ef['T'] + basecount_ef['G'])
	reverse = basecount_nc['r'] + basecount_nc['a'] + basecount_nc['c'] + basecount_nc['t'] + basecount_nc['g'] + (basecount_ef['R'] + basecount_ef['A'] + basecount_ef['C'] + basecount_ef['T'] + basecount_ef['G'] + basecount_ef['r'] + basecount_ef['a'] + basecount_ef['c'] + basecount_ef['t'] + basecount_ef['g'])
	forward_var = basecount_nc[var] + basecount_ef[var.lower()] + basecount_ef[var]
	reverse_var = basecount_nc[var.lower()] + basecount_ef[var] + basecount_ef[var.lower()]
	if reverse_var == 0 and forward_var == 0:
		return "F"
	if reverse_var == 0 and reverse == 0:
		return "T"
	elif reverse_var == 0 and reverse != 0:
		var_ratio = float(reverse_var)/float(forward_var)
		ratio = float(reverse)/float(forward)
#		print "merge reverse var 0", binom.pmf(forward_var, forward_var + reverse_var, float(forward)/(float(forward)+float(reverse))), reverse, forward, reverse_var, forward_var
		if binom.pmf(forward_var, forward_var + reverse_var, float(forward)/(float(forward)+float(reverse))) < STRAND_BIAS_BINOMIAL_PROB:
			return "F"
		else:
			return "T" 
#		if var_ratio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL:
#			return "F"
#		else:
#			return "T"
	elif forward_var == 0 and forward == 0:
		return "T"
	elif forward_var == 0 and forward != 0:
		var_ratio = float(forward_var)/float(reverse_var)
		ratio = float(forward)/float(reverse)
#		print "merge forward var 0", binom.pmf(reverse_var, forward_var + reverse_var, float(reverse)/(float(forward)+float(reverse))), reverse, forward, reverse_var, forward_var
		if binom.pmf(reverse_var, forward_var + reverse_var, float(reverse)/(float(forward)+float(reverse))) < STRAND_BIAS_BINOMIAL_PROB:
			return "F"
		else:
			return "T"
#		if var_ratio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL:
#			return "F"
#		else:
#			return "T"
	else:
		var_ratio = float(forward_var)/float(reverse_var)
		ratio = float(forward)/float(reverse)
		if var_ratio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL or ratio/var_ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL:
			return "F"
		else:
			return "T"



def filter_both_strand_above_average_merge(basecount_notCombined, basecount_extendedFrags):
# input: base counts, non-reference alleles
# output: if the locus passes two filters
	# The first filter: the variant allele should be observed on both strands
	# The second filter: the proportion of the variant alelle among all non-reference alleles should above a certain threshold of average 
# related global variables:
	# THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF
	if_both_observed = 0
	if_above_average = 0
	for base in ['A', 'C', 'T', 'G']:
		if ( basecount_notCombined[base.lower()] > 0  and basecount_notCombined[base] > 0 ) or basecount_extendedFrags[base.lower()] > 0 or basecount_extendedFrags[base] > 0 :
			if_both_observed = 1
		if ((basecount_notCombined[base] + basecount_notCombined[base.lower()]) + 2*basecount_extendedFrags[base.lower()] + 2*basecount_extendedFrags[base] ) >= THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF*((sum(basecount_notCombined.values()) + 2*sum(basecount_extendedFrags.values()) - basecount_notCombined['R'] - basecount_notCombined['r'] - 2*basecount_extendedFrags['R'] - 2*basecount_extendedFrags['r'])):
			if_above_average = 1
	if if_both_observed == 1:
		both_observed = "T"
	else:
		both_observed = "F"
	if if_above_average == 1:
		above_average = "T"
	else:
		above_average = "F"
	return both_observed, above_average




def filter_strand_bias_unmerge(basecount, var):
# input: base counts
# output: if the locus passes the filter
		# the percentage of reads from forward strand
# related global variables:
	# TOTAL_STRAND_BIAS_RATIO_HIGH
	# TOTAL_STRAND_BIAS_RATIO_LOW
	# print basecount_nc, basecount_ef
	forward = basecount['R'] + basecount['A'] + basecount['C'] + basecount['T'] + basecount['G'] 
	reverse = basecount['r'] + basecount['a'] + basecount['c'] + basecount['t'] + basecount['g']
	forward_var = basecount[var] 
	reverse_var = basecount[var.lower()]
	if reverse_var == 0 and reverse == 0:
		return "T"
	elif reverse_var == 0 and reverse != 0:
		var_ratio = float(reverse_var)/float(forward_var)
		ratio = float(reverse)/float(forward)
#		print "unmerge reverse var 0", binom.pmf(forward_var, forward_var + reverse_var, float(forward)/(float(forward)+float(reverse))), reverse, forward, reverse_var, forward_var
		if binom.pmf(forward_var, forward_var + reverse_var, float(forward)/(float(forward)+float(reverse))) < STRAND_BIAS_BINOMIAL_PROB:
			return "F"
		else:
			return "T"
		if binom.cdf(reverse_var, forward_var + reverse_var, float(reverse)/(float(forward)+float(reverse))) < STRAND_BIAS_BINOMIAL_PROB:
			return "F"
		else:
			return "T"
#		if var_ratio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL:
#			return "F"
#		else:
#			return "T"
	elif forward_var == 0 and forward == 0:
		return "T"
	elif forward_var == 0 and forward != 0:
		var_ratio = float(forward_var)/float(reverse_var)
		ratio = float(forward)/float(reverse)
#		print "unmerge forward var 0", binom.pmf(reverse_var, forward_var + reverse_var, float(reverse)/(float(forward)+float(reverse))), reverse, forward, reverse_var, forward_var
		if binom.pmf(reverse_var, forward_var + reverse_var, float(reverse)/(float(forward)+float(reverse))) < STRAND_BIAS_BINOMIAL_PROB:
			return "F"
		else:
			return "T"
		if binom.cdf(forward_var, forward_var + reverse_var, float(forward)/(float(forward)+float(reverse))) < STRAND_BIAS_BINOMIAL_PROB:
			return "F"
		else:
			return "T"
#		if var_ratio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL:
#			return "F"
#		else:
#			return "T"
	else:
		var_ratio = float(forward_var)/float(reverse_var)
		ratio = float(forward)/float(reverse)
	#	print forward_var, reverse_var, forward, reverse, var_ratio, ratio
		if var_ratio/ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL or ratio/var_ratio > STRAND_BIAS_RATIO_VARIANT_TO_ALL:
			return "F"
		else:
			return "T"



def filter_both_strand_above_average_unmerge(basecount):
# input: base counts, non-reference alleles
# output: if the locus passes two filters
	# The first filter: the variant allele should be observed on both strands
	# The second filter: the proportion of the variant alelle among all non-reference alleles should above a certain threshold of average 
# related global variables:
	# THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF
	if_both_observed = 0
	if_above_average = 0
	for base in ['A', 'C', 'T', 'G']:
		if basecount[base.lower()] > 0  and basecount[base] > 0 :
			if_both_observed = 1
		if (basecount[base] + basecount[base.lower()]) >= THRESHOLD_VARIANT_ALLELE_PROPORTION_IN_NONREF*(sum(basecount.values()) - basecount['R'] - basecount['r']):
			if_above_average = 1
	if if_both_observed == 1:
		both_observed = "T"
	else:
		both_observed = "F"
	if if_above_average == 1:
		above_average = "T"
	else:
		above_average = "F"
	return both_observed, above_average



def filter_all_mapping_quality(maplist):
	if sum(maplist < OK_MAPQUAL) > THRESHOLD_OK_MAPQUAL_COUNT_ALL:
		is_all_mapping_qual = 'T'
	elif (sum(maplist > LOW_MAPQUAL) > LOW_MAPQUAL_COUNT_FOR_DETECTION_ALL and len(maplist) < DEPTH_ACTIVATE_LOW_MAPQUAL_COUNT_ALL) or sum(maplist > LOW_MAPQUAL)/len(maplist) > LOW_MAPQUAL_FRAC_FOR_DETECTION:
		is_all_mapping_qual = 'F'
	elif np.median(maplist) > AVG_MAPQUAL_FOR_DETECTION:
		is_all_mapping_qual = 'F'
	else:
		is_all_mapping_qual = 'T'
	return is_all_mapping_qual



def filter_variant_mapping_quality(maplist, basestring, variant_base):
	baselist = np.array(list(basestring.upper()))
	var_id = (baselist == variant_base)
	variant_maplist = maplist[baselist == variant_base]
	if len(variant_maplist) == 0:
		is_variant_mapping_qual = 'F'
	elif sum(variant_maplist <= OK_MAPQUAL) > THRESHOLD_OK_MAPQUAL_COUNT:
		is_variant_mapping_qual = 'T'
	elif sum(variant_maplist > LOW_MAPQUAL) > LOW_MAPQUAL_COUNT_FOR_DETECTION_VARIANT and len(variant_maplist) < DEPTH_ACTIVATE_LOW_MAPQUAL_COUNT:
		is_variant_mapping_qual = 'F'
	elif sum(variant_maplist > LOW_MAPQUAL)/len(variant_maplist) > LOW_MAPQUAL_FRAC_FOR_DETECTION:
		is_variant_mapping_qual = 'F'
	elif np.median(variant_maplist) > AVG_MAPQUAL_FOR_DETECTION:
		is_variant_mapping_qual = 'F'
	else:
		is_variant_mapping_qual = 'T'
	return is_variant_mapping_qual



def filter_variant_base_quality(quallist, basestring, variant_base):
	baselist = np.array(list(basestring.upper()))
	var_id = (baselist == variant_base)
	ref_id = (baselist == 'R')
	nonref_id = (baselist != variant_base) * (baselist != 'R')
	if sum(var_id) == 0:
		return 'F'
	if sum(ref_id) == 0:
		return 'T'
	avg_var_qual = np.median(quallist[var_id])
	avg_ref_qual = np.median(quallist[ref_id])
	if avg_var_qual/avg_ref_qual > THRESHOLD_VARIANT_QUALITY_TO_REFERENCE_QUALITY and sum(quallist[var_id] < OK_BASEQUAL) < THRESHOLD_OK_BASEQUAL_COUNT:
		return 'F'
	else:
		return 'T'



#dbSNP_file = "/u/project/xjzhou/shuoli/variants_public_database/dbSNP/common_all_20170710_shuo_snvonly.txt"
dbSNP = defaultdict(list)

def import_dbSNP(dbSNP_file):
	f = open(dbSNP_file, 'r')
	for line in f:
		if line[0] == "#":
			continue
		sp = line.strip().split('\t')
		if "VC=SNV" not in sp[7]:
			continue
		spp = sp[2].split(',')
		for i in spp:
			if len(i) > 1:
				continue
			if i not in "ACGT":
				dbSNP["chr" + sp[0] + '-' + 'A'].append(sp[1])
				dbSNP["chr" + sp[0] + '-' + 'C'].append(sp[1])
				dbSNP["chr" + sp[0] + '-' + 'G'].append(sp[1])
				dbSNP["chr" + sp[0] + '-' + 'T'].append(sp[1])
			else:
				dbSNP["chr" + sp[0] + '-' + i].append(sp[1])
	f.close()
	return


def filter_dbSNP(chrom, position, variant_base):
	if position in dbSNP[chrom + '-' + variant_base]:
		return 'F'
	else:
		return 'T'



def filter_supporting_reads_unmerge(basestring, variant_base):
	if basestring.count(variant_base) + basestring.count(variant_base.lower()) < COUNT_FOR_STRONG_EVIDENCE:
		return 'F'
	else:
		return 'T'



def filter_supporting_fragment_merge(basestring_merge, variant_base):
	if basestring_merge.count(variant_base) + basestring_merge.count(variant_base.lower()) < COUNT_FOR_STRONG_EVIDENCE:
		return 'F'
	else:
		return 'T'


def filter_VAF(VAF, depth):
	VAF_FOR_STRONG_EVIDENCE = 3.0/depth
	if VAF < VAF_FOR_STRONG_EVIDENCE:
		return 'F', "weak VAF"
	else:
		return 'T', "strong VAF"


def final_estimation_with_Jenks(VAF_list):
	if len(VAF_list) == 0:
		return 0.0, 0, 0
	if len(VAF_list) < 20:
		return sum(VAF_list)/float(len(VAF_list)), len(VAF_list), 0
	tmp1 = 0
	max_iter_lim = 10
	min_iter_lim = 3
	max_iter_x_lim = int(round(float(len(VAF_list))/20.0)) ###?
	max_lim = max(min_iter_lim, min(max_iter_lim, max_iter_x_lim)) + 1
	for ct in range(min_iter_lim, max_lim):
		groups, breaks = jenks(VAF_list, ct)
		tmp2 = float(sum(groups[len(groups)-1]))/len(groups[len(groups)-1])
		if tmp1 == tmp2:
			break
		else:
			tmp1 = tmp2
	include = []
	n = 0
	for i in range(len(groups)):
		group = groups[len(groups)-i-1]
		include += group
		n += 1
		if len(include) >= VARIANT_COUNT_STOP_ADDING_IN_JENK_ESTIMATION:
			break
	if len(include) > VARIANT_COUNT_FOR_JENK_ESTIMATION:
		include = include[0:VARIANT_COUNT_FOR_JENK_ESTIMATION]
	return float(sum(include))/len(include), len(include), n


def filter_supporting_count_by_binomial_test(jenks_estimate, basestring, variant_base):
	cnt = len(basestring)
	var_cnt = basestring.upper().count(variant_base)
	p = binom.cdf(var_cnt, cnt, jenks_estimate/2)
	if p < P_VALUE_BINOMIAL_TEST_JENKS_ESTIMATE:
		return 'F'
	else:
		return 'T'



def filter_normal_coverage(count_normal):
	if count_normal <= 7:
		return 'F' # low confidence in normal
	if count_normal > 7:
		return 'T'


def filter_tumor_coverage(count):
	if count <= 10:
		return 'F', "low confidence call in plasma"
	if count > 10 and count <= 50:
		return 'T', "low confidence VAF in plasma"
	if count > 50:
		return 'T', "high confidence in plasma"




