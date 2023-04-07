from parameter import *
from _probability import *
from _filter import *
import numpy as np
import sys


# POSTERIOR PROBABILITY CALCULATION

# log( P( G, X | theta ) )
# P( G | X , theta ) ~ P( G, X | theta )
def calculate_joint_genotype_logposterior(tumor_fraction, basestring_merge, quallist_merge, maplist_merge, basestring_normal, quallist_normal, maplist_normal, variant_base):
# input: tumor fraction, converted basestring, base quality list
# output: log posterior of all joint genotypes log( P( G, X | theta ) ), non-reference alleles, base count string
# P( G | X , theta ) is computed, but P( G | X , theta ) ~ P( G, X | theta )
	joint_genotype_list = list(PRIOR.keys())
	logposterior_dict = {}
	for i in range(len(joint_genotype_list)):
		#logposterior_dict[joint_genotype_list[i]] = calculate_joint_genotype_tumor_fraction_loglikelihood(tumor_fraction, basestring_merge, quallist_merge, maplist_merge, joint_genotype_list[i], variant_base) + calculate_joint_genotype_tumor_fraction_loglikelihood(0.0, basestring_normal, quallist_normal, maplist_normal, joint_genotype_list[i], variant_base)  #  + math.log(PRIOR[joint_genotype_list[i]])) ####
		logposterior_dict[joint_genotype_list[i]] = (len(basestring_merge) + len(basestring_normal)) * (1.0/len(basestring_merge) * calculate_joint_genotype_tumor_fraction_loglikelihood(tumor_fraction, basestring_merge, quallist_merge, maplist_merge, joint_genotype_list[i], variant_base) + 1.0/len(basestring_normal) * calculate_joint_genotype_tumor_fraction_loglikelihood(0.0, basestring_normal, quallist_normal, maplist_normal, joint_genotype_list[i], variant_base))  #  + math.log(PRIOR[joint_genotype_list[i]])) ####
	return logposterior_dict

def calculate_confidence(logposterior_dict):
# input: log posterior list of all joint genotypes
# output: confidence score of the genotype with the maximum posterior
# confidence score is computed by the ratio of the maximum posterior and the second maximum posterior
	max_var = max(logposterior_dict["AA/AB"], logposterior_dict["AA/BB"])
	max_alt = max([logposterior_dict[i] for i in PRIOR.keys() if i != "AA/AB" and i != "AA/BB"])
	confidence = convert_to_Decimal(max_var).exp()/convert_to_Decimal(max_alt).exp()
	return confidence

def genotype(logposterior_dict):
# input: log posterior list of all joint genotypes, non-reference alleles
# output: joint genotype, non-reference alleles
# joint genotype is determined by maximizing a posterior
	return max(logposterior_dict, key=logposterior_dict.get)


# REPORT
def output_genotyping(output, chrom, pos, ref, cnt_normal, cnt_var_normal, cnt_tumor, cnt_var_tumor, cnt_tumor_merge, cnt_var_tumor_merge, variant_base, joint_genotype_MAP, VAF_tumor, confidence, joint_genotype_MAP_merge, VAF_tumor_merge, confidence_merge, baseinfo_list):
# input: chromosome, position, basestring, joint genotype, non-reference alleles, base count string, confidence, output file
# output: none
# only somatic variants are output currently
# all variants should be output
# related global variables:
	# DEPTH_FOR_DETECTION
	res = [chrom, pos, ref.upper(), cnt_normal, cnt_var_normal, cnt_tumor, cnt_var_tumor, cnt_tumor_merge, cnt_var_tumor_merge, variant_base, joint_genotype_MAP, VAF_tumor, confidence, joint_genotype_MAP_merge, VAF_tumor_merge, confidence_merge] + baseinfo_list
	res = ["\"" + str(i) + "\"" for i in res]
	output.write('\t'.join(res) + '\n')
	return

def call_variants(finput, foutput_call, tumor_fraction):
# input: input filename, variant calling output filename, tumor fraction
# output: none
	# unfiltered variants are stored in foutput_call
	f = open(finput, 'r')
	n = 0
	m = 0
	output = open(foutput_call, 'w')
	for line in f:
		m += 1
		sp = line.split('\t')
		sp[len(sp)-1] = sp[len(sp)-1].strip()
		chrom = sp[0]
		pos = int(sp[1])
		ref = sp[2]
		nread = len(sp[19])
		basestring = sp[19]
		basestring_all = sp[4]
		if len(basestring) == 0:
			continue
		quallist = string_to_qual(sp[20])
		maplist = string_to_qual(sp[21])
		basestring_normal = sp[22]
		quallist_normal = string_to_qual(sp[23])
		maplist_normal = string_to_qual(sp[24])
		basestring_extendedFrags = sp[25]
		quallist_extendedFrags = string_to_qual(sp[26])
		maplist_extendedFrags = string_to_qual(sp[27])
		if len(sp) <= 29:
			print(line)
		basestring_notCombined = sp[28]
		quallist_notCombined = string_to_qual(sp[29])
		maplist_notCombined = string_to_qual(sp[30])
		basestring_normal_all = sp[8]
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
		if basestring_normal_all.upper().count(variant_base_all) > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER:
			continue
		if basestring_normal_all.upper().count(variant_base) > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER:
			continue
		if basestring_normal.upper().count(variant_base) > GERMLINE_VARIANT_COUNT:
			continue
		basestring_merge = basestring_extendedFrags + basestring_notCombined
		quallist_merge = np.concatenate([quallist_extendedFrags, quallist_notCombined])
		maplist_merge = np.concatenate([maplist_extendedFrags, maplist_notCombined])
		if len(basestring_normal) < DEPTH_FOR_DETECTION_NORMAL:
			continue
		if len(basestring_merge) < DEPTH_FOR_DETECTION_TUMOR or len(basestring) < DEPTH_FOR_DETECTION_TUMOR:
			continue
		cnt_var_tumor = basestring.upper().count(variant_base)
		cnt_var_tumor_merge = basestring_merge.upper().count(variant_base)
		cnt_var_normal = basestring_normal.upper().count(variant_base)
		VAF_tumor_merge = float(cnt_var_tumor_merge)/float(len(basestring_merge))
		VAF_tumor = float(cnt_var_tumor)/float(len(basestring))
		VAF_normal = float(cnt_var_normal)/float(len(basestring_normal))
		cnt_alt_normal = sum([1 for i in basestring_normal.upper() if i != 'R'])
		cnt_normal = len(basestring_normal)
		#print VAF_normal, cnt_var_normal
		if VAF_normal > SOMATIC_VAF_THRESHOLD_IN_NORMAL:
			continue
		if cnt_var_normal > 0 and binom.cdf(cnt_var_normal, len(basestring_normal), VAF_tumor) > NORMAL_COUNT_BINOM:
			continue	
		if cnt_var_normal > NORMAL_COUNT_VAR or cnt_alt_normal > NORMAL_COUNT_ALT:
			continue
		if cnt_var_normal > 0 and (VAF_normal >= VAF_tumor_merge or VAF_normal >= VAF_tumor or cnt_var_normal >= cnt_var_tumor or cnt_var_normal >= cnt_var_tumor_merge):
			continue
		#print chrom, pos
		logposterior_dict_merge = calculate_joint_genotype_logposterior(tumor_fraction, basestring_merge, quallist_merge, maplist_merge, basestring_normal, quallist_normal, maplist_normal, variant_base)
		confidence_merge = calculate_confidence(logposterior_dict_merge)
		joint_genotype_MAP_merge = genotype(logposterior_dict_merge)
		logposterior_dict = calculate_joint_genotype_logposterior(tumor_fraction, basestring, quallist, maplist, basestring_normal, quallist_normal, maplist_normal, variant_base)
		confidence = calculate_confidence(logposterior_dict)
		joint_genotype_MAP = genotype(logposterior_dict)
		#print joint_genotype_MAP
		if joint_genotype_MAP == "AA/BB" or joint_genotype_MAP == "AA/AB" or joint_genotype_MAP_merge == "AA/BB" or joint_genotype_MAP_merge == "AA/AB":
			n += 1
			if n % 100 == 0:
				print(n, m)
			if VAF_tumor <= GENOTYPE_THRESHOLD and VAF_tumor_merge <= GENOTYPE_THRESHOLD:
				output_genotyping(output, chrom, pos, ref, cnt_normal, cnt_var_normal, len(basestring), cnt_var_tumor, len(basestring_merge), cnt_var_tumor_merge, variant_base, joint_genotype_MAP, VAF_tumor, confidence, joint_genotype_MAP_merge, VAF_tumor_merge, confidence_merge, sp[19:31])
	f.close()
	output.close()


if __name__ == "__main__":
	finput = sys.argv[1]
	foutput_call = sys.argv[2]
	tumor_fraction = float(sys.argv[3])
	MERGED_VAF_THRESHOLD = float(sys.argv[4])
	depth = float(sys.argv[5])
	GENOTYPE_THRESHOLD = MERGED_VAF_THRESHOLD * 1.2
	GERMLINE_VARIANT_COUNT = get_germline_variant_count_threshold(depth)
	NORMAL_COUNT_BINOM = get_normal_count_binom_threshold(depth)
	call_variants(finput, foutput_call, tumor_fraction)

