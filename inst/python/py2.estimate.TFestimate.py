from parameter import *
from _filter import *
from _probability import *
import numpy as np
import sys

# PREPROCESSING WITH PREDEFINED HOTSPOT LIST
def import_predefined_hotspot_prior(filename):
# input file format:	VCF or BED (chromosome	position) 
					# tab deliminated
# related global variables:
	# HOTSPOT_PRIOR
	f = open(filename, 'r')
	for line in f:
		sp = line.strip().split('\t')
		HOTSPOT_PRIOR.append((sp[0], int(sp[1])))
	return


def select_hotspot(chrom, pos, basestring, quallist, maplist, basestring_normal, quallist_normal, maplist_normal,  basestring_extendedFrags, quallist_extendedFrags, maplist_extendedFrags, basestring_notCombined, quallist_notCombined, maplist_notCombined, basecount_notCombined, basecount_extendedFrags, variant_base):
# input: chrom, position, number of reads covering the position, 
 		  # basestring after rebase, quallist after rebase and convert
# output: none (add position to HOTSPOT if meet requirement)
# determine if a locus can be involved in hotspot set
# related global parameters:
	# HOTSPOT_VARIANT_QUALITY_MAX_PERCENTAGE
	# GERMLINE_VAF_THRESHOLD
	# HOTSPOT
	# require DEPTH_FOR_ESTIMATION reads
	if nread < DEPTH_FOR_ESTIMATION:
		return
	basecount = count_base(basestring)
	basestring_merge = basestring_extendedFrags + basestring_notCombined
	quallist_merge = np.concatenate([quallist_extendedFrags, quallist_notCombined])
	maplist_merge = np.concatenate([maplist_extendedFrags, maplist_notCombined])
	if basestring_normal.upper().count(variant_base) > GERMLINE_VARIANT_COUNT:
		return
	# require DEPTH_FOR_ESTIMATION fragments
	if len(basestring_merge) < DEPTH_FOR_ESTIMATION:
		return
	# exclude loci without enough germline information
	if len(basestring_normal) < DEPTH_FOR_GERMLINE:
		return
	# require high quality variant alleles
	if len(basestring) == 0 or len(basestring_merge) == 0 or len(basestring_normal) == 0:
		return
	cnt_non_ref_read = sum([1 for i in basestring.upper() if i == variant_base])
	mean_qual = float(sum(quallist))/float(len(quallist))
	if float(cnt_non_ref_read) < mean_qual*float(len(basestring)):
		return	
	cnt_non_ref_read_merge = sum([1 for i in basestring_merge if i.upper() == variant_base])
	mean_qual_merge = np.mean(quallist_merge)
	if float(cnt_non_ref_read_merge) < mean_qual_merge*float(len(basestring_merge)):
		return	
	# merged VAF lower than threshold
	cnt_var_read_merge = sum([1 for i in basestring_merge if i.upper() == variant_base])
	if float(cnt_var_read_merge)/float(len(basestring_merge)) > MERGED_VAF_THRESHOLD:
		return
	# exclude possible germline variants
	cnt_var_normal_read = sum([1 for i in basestring_normal if i.upper() == variant_base])
	if float(cnt_var_normal_read)/float(len(basestring_normal)) > GERMLINE_VAF_THRESHOLD_IN_NORMAL:
		return
	# require high mapping quality in both tumor and normal
	if sum(maplist == 0) > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION or sum(maplist_normal == 0) > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION or sum(maplist_merge == 0) > ZERO_MAPQUAL_COUNT_FOR_ESTIMATION:
		return
	if np.mean(maplist[maplist > 0]) > AVG_MAPQUAL_FOR_ESTIMATION or np.mean(maplist_normal[maplist_normal > 0]) > AVG_MAPQUAL_FOR_ESTIMATION or np.mean(maplist_merge[maplist_merge > 0]) > AVG_MAPQUAL_FOR_ESTIMATION:
		return
	# exclude strand bias
	strand_bias_merge = filter_strand_bias_merge(basecount_notCombined, basecount_extendedFrags, variant_base)
	# require enough percentage in non-reference reads
	both_observed_merge, above_average_merge = filter_both_strand_above_average_merge(basecount_notCombined, basecount_extendedFrags)
	supporting_frag_count_merge = filter_supporting_fragment_merge(basestring_merge, variant_base)
	strand_bias_unmerge = filter_strand_bias_unmerge(basecount, variant_base)
	# require enough percentage in non-reference reads
	both_observed_unmerge, above_average_unmerge = filter_both_strand_above_average_unmerge(basecount)
	supporting_read_count_unmerge = filter_supporting_reads_unmerge(basestring, variant_base)
	if strand_bias_merge != "T" or both_observed_merge != "T" or above_average_merge != "T" or supporting_frag_count_merge != "T":
		return
	if strand_bias_unmerge != "T" or both_observed_unmerge != "T" or above_average_unmerge != "T" or supporting_read_count_unmerge != "T" :
		return
	HOTSPOT.append((mean_qual_merge, cnt_var_normal_read, float(cnt_var_read_merge)/float(len(basestring_merge)), chrom, pos, basestring, quallist, maplist, basestring_normal, quallist_normal, maplist_normal,  basestring_extendedFrags, quallist_extendedFrags, maplist_extendedFrags, basestring_notCombined, quallist_notCombined, maplist_notCombined, variant_base))
	return


def select_hotspot_from_predifined_hotspot_prior(chrom, pos, basestring, quallist, maplist, basestring_normal, quallist_normal, maplist_normal,  basestring_extendedFrags, quallist_extendedFrags, maplist_extendedFrags, basestring_notCombined, quallist_notCombined, maplist_notCombined, basecount_notCombined, basecount_extendedFrags, variant_base):
# input: chrom, position, number of reads covering the position, 
 		  # basestring after rebase, quallist after rebase and convert
# output: none (add position to HOTSPOT if meet requirement)
# determine if a locus in hotspot prior can be involved in hotspot set
# call select_hotspot
	if (chrom, pos) in HOTSPOT_PRIOR:
		print (chrom, pos)
		select_hotspot(chrom, pos, basestring, quallist, maplist, basestring_normal, quallist_normal, maplist_normal,  basestring_extendedFrags, quallist_extendedFrags, maplist_extendedFrags, basestring_notCombined, quallist_notCombined, maplist_notCombined, basecount_notCombined, basecount_extendedFrags, variant_base)
	return



def calculate_tumor_fraction_likelihood(tumor_fraction, basestring_merge, quallist_merge, maplist_merge, basestring_normal, quallist_normal, maplist_normal, variant_base):
# input: tumor fraction, converted basestring, base quality list
# output: likelihood P( X | theta )
	joint_genotype_list = PRIOR.keys()
	loglikelihood_dict = dict([])
	likelihood = 0
	for i in joint_genotype_list:
		loglikelihood_dict[i] = calculate_joint_genotype_tumor_fraction_loglikelihood(tumor_fraction, basestring_merge, quallist_merge, maplist_merge, i, variant_base) + calculate_joint_genotype_tumor_fraction_loglikelihood(0.0, basestring_normal, quallist_normal, maplist_normal, i, variant_base)
		likelihood += convert_to_Decimal((loglikelihood_dict[i]+math.log(EST_PRIOR[i]))).exp()
	return likelihood


def estimate_tumor_fraction(HOTSPOT):
# with pileup file (current version)
# input: all HOTSPOTS
# output: estimated tumor fraction
# tumor fraction is estimated with maximum likelihood
	ratio = [GRIDWIDTH * i for i in range(MAXSEARCH_WITH_NORMAL)]
	tumor_fraction_likelihood = [convert_to_Decimal(0) for i in range(MAXSEARCH_WITH_NORMAL)]
	for ratioind in range(MAXSEARCH_WITH_NORMAL):
	#	print ratioind
		for spot in HOTSPOT:
			tumor_fraction_likelihood[ratioind] += (calculate_tumor_fraction_likelihood(ratio[ratioind], spot[11]+spot[14], np.concatenate([spot[12], spot[15]]), np.concatenate([spot[13], spot[16]]), spot[8], spot[9], spot[10], spot[17]).ln())/(len(spot[11]) + len(spot[14]) + len(spot[8]))
	estind = tumor_fraction_likelihood.index(max(tumor_fraction_likelihood))
	est = ratio[estind]
	return est, tumor_fraction_likelihood


if __name__ == "__main__":
	filename = sys.argv[1]
	MERGED_VAF_THRESHOLD = float(sys.argv[2])
	file_prefix = sys.argv[3]
	depth = float(sys.argv[4])
	VAF_output = sys.argv[5]
	estimate_output = sys.argv[6]
	GERMLINE_VARIANT_COUNT = get_germline_variant_count_threshold(depth)
	f = open(filename, 'r')
	for line in f:
		sp = line.strip().split('\t')
		chrom = sp[0]
		if "X" in chrom or "M" in chrom or "Y" in chrom:
			continue
		pos = sp[1]
		nread = len(sp[19]) 
		basestring_all = sp[4]
		basestring = sp[19]
		if len(basestring) == 0:
			continue
		quallist = string_to_qual(sp[20])
		maplist = string_to_qual(sp[21])
		if len(sp) < 31:
			continue
		#print line 
		basestring_normal_all = sp[8]
		basestring_normal = sp[22]
		quallist_normal = string_to_qual(sp[23])
		maplist_normal = string_to_qual(sp[24])
		basestring_extendedFrags = sp[25]
		quallist_extendedFrags = string_to_qual(sp[26])
		maplist_extendedFrags = string_to_qual(sp[27])
		basestring_notCombined = sp[28]
		quallist_notCombined = string_to_qual(sp[29])
		maplist_notCombined = string_to_qual(sp[30])
		basecount_notCombined_all = count_base(sp[12])
		basecount_extendedFrags_all = count_base(sp[16])
		variant_base_all = find_major_variant(basecount_notCombined_all, basecount_extendedFrags_all)
		basecount_notCombined = count_base(basestring_notCombined)
		basecount_extendedFrags = count_base(basestring_extendedFrags)
		variant_base = find_major_variant(basecount_notCombined, basecount_extendedFrags)
		basecount_tumor_all = count_base(basestring_all)
		basecount_tumor = count_base(basestring)
		triallele_tumor_all = filter_triallelic_position(basecount_tumor_all, variant_base_all, depth)
		triallele_tumor = filter_triallelic_position(basecount_tumor, variant_base, depth)
		if variant_base != variant_base_all:
			continue
		if basestring.count(variant_base.upper()) + basestring.count(variant_base.lower()) == 0:
			continue
		if triallele_tumor_all == "F" or triallele_tumor == "F":
			continue
		if basestring_normal_all.upper().count(variant_base) > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER:
			continue
		if basestring_normal_all.upper().count(variant_base_all) > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER:
			continue
		select_hotspot(chrom, pos, basestring, quallist, maplist, basestring_normal, quallist_normal, maplist_normal,  basestring_extendedFrags, quallist_extendedFrags, maplist_extendedFrags, basestring_notCombined, quallist_notCombined, maplist_notCombined, basecount_notCombined, basecount_extendedFrags, variant_base)
	f.close()
	if len(HOTSPOT) > MAX_NUMBER_OF_HOTSPOT:
		from operator import itemgetter, attrgetter
		a = sorted(HOTSPOT, key=itemgetter(0))
		a = sorted(a, key=itemgetter(2), reverse = True)
		a = sorted(a, key=itemgetter(1))
		#HOTSPOT_ALL = [i for i in HOTSPOT]
		HOTSPOT = a[0:MAX_NUMBER_OF_HOTSPOT] #####
	est, tumor_fraction_likelihood = estimate_tumor_fraction(HOTSPOT)
#	print est
	VAF = []
	for i in HOTSPOT:
		vaf = sum([1.0 for j in i[11]+i[14] if j.upper() == i[17]])/float(len(i[11]+i[14]))
		VAF.append(vaf)
	VAF.sort()
	VAF
	h = open(VAF_output, 'w')
	for i in VAF:
		h.write(str(i) + '\n')

	h.close()
#	est, tumor_fraction_likelihood = estimate_tumor_fraction(HOTSPOT)
#	print est
	g = open(estimate_output, 'w')
	g.write(str(est) + '\n')
	g.write("====================" + '\n')
	for i in tumor_fraction_likelihood:
		g.write(str(i) + '\n')
	g.close()



