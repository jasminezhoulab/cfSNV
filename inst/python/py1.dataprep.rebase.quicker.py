import sys
from parameter import *
from _probability import *
from _filter import *
import numpy as np


# PREPROCESSING WITH PILEUP FILE

def rebase_with_frequent_variant_allele(string, qual_string, map_string):
# input: a line of pileup file
# output: list of [base, qual] pairs 
# 1. find most frequent variant allele V
# 2. extract bases related to V and ref allele
# 3. extract quals corresponding to bases
# 4. multiple most frequent variant alleles are allowed
# 5. consistency check of #bases and #quals selected in the end 
	mod_base_string = ''
	mod_qual_string = ''
	mod_map_string = ''
	mod_base_string_pass = ''
	mod_qual_string_pass = ''
	mod_map_string_pass = ''
	cdict = dict([])
	insert = 0
	delete = 0
	cdict['A'] = ['', '', '', '', '', '']
	cdict['C'] = ['', '', '', '', '', '']
	cdict['G'] = ['', '', '', '', '', '']
	cdict['T'] = ['', '', '', '', '', '']
	i = 0
	j = 0
	# convert string and qual_string
	while i < len(string):
   #     print i
		base = string[i]
		if base in '.': #match on forward
			mod_base_string += 'R'
			mod_qual_string += qual_string[j]
			mod_map_string += map_string[j]
			j += 1
			i += 1
			if 10.0**(-float(ord(qual_string[j-1])-33)/10.0) > BASEQUAL_THRESHOLD:
				continue
			if 10.0**(-float(ord(map_string[j-1])-33)/10.0) > MAPQUAL_THRESHOLD:
				continue
			mod_base_string_pass += 'R'
			mod_qual_string_pass += qual_string[j-1]
			mod_map_string_pass += map_string[j-1]
			continue
		if base in ',': #match on reverse
			mod_base_string += 'r'
			mod_qual_string += qual_string[j]
			mod_map_string += map_string[j]
			j += 1
			i += 1
			if 10.0**(-float(ord(qual_string[j-1])-33)/10.0) > BASEQUAL_THRESHOLD:
				continue
			if 10.0**(-float(ord(map_string[j-1])-33)/10.0) > MAPQUAL_THRESHOLD:
				continue
			mod_base_string_pass += 'r'
			mod_qual_string_pass += qual_string[j-1]
			mod_map_string_pass += map_string[j-1]
			continue
		if base in 'ACGT': #mismatch on forward
			cdict[base.upper()][0] += base
			cdict[base.upper()][1] += qual_string[j]
			cdict[base.upper()][2] += map_string[j]
			j += 1
			i += 1
			if 10.0**(-float(ord(qual_string[j-1])-33)/10.0) > BASEQUAL_THRESHOLD:
				continue
			if 10.0**(-float(ord(map_string[j-1])-33)/10.0) > MAPQUAL_THRESHOLD:
				continue
			cdict[base.upper()][3] += base
			cdict[base.upper()][4] += qual_string[j-1]
			cdict[base.upper()][5] += map_string[j-1]
			continue
		if base in 'acgt': #mismatch on reverse
			cdict[base.upper()][0] += base
			cdict[base.upper()][1] += qual_string[j]
			cdict[base.upper()][2] += map_string[j]
			j += 1
			i += 1
			if 10.0**(-float(ord(qual_string[j-1])-33)/10.0) > BASEQUAL_THRESHOLD:
				continue
			if 10.0**(-float(ord(map_string[j-1])-33)/10.0) > MAPQUAL_THRESHOLD:
				continue
			cdict[base.upper()][3] += base
			cdict[base.upper()][4] += qual_string[j-1]
			cdict[base.upper()][5] += map_string[j-1]
			continue
		if base in 'Nn': #mismatch on forward
			j += 1
			i += 1
			continue
		if base in '*':
			i += 1
			j += 1
			delete += 1
		if base == '^': #start of a read + a phred score
			i += 2
			continue
		if base == '$': #end of a read + a phred score
			i += 1
			continue
		if base in '+-': #indel
			if base == '+':
				insert += 1
			else:
				delete += 1
			i += 1
			base = string[i]
			num = ''
			while base in '0123456789':
				num += base
				i += 1
				base = string[i]
			num = int(num) #num of indel bases
			i += num
			continue	
	# find max frequent variant allele
	length = [len(cdict['A'][0]), len(cdict['C'][0]), len(cdict['G'][0]), len(cdict['T'][0])]
	max_length = max(length)
	if max_length == 0:
		return ([mod_base_string, mod_qual_string, mod_map_string, mod_base_string_pass, mod_qual_string_pass, mod_map_string_pass], insert, delete, len(mod_base_string) + insert + delete)
	new_base_string = mod_base_string
	new_qual_string = mod_qual_string
	new_map_string = mod_map_string
	new_base_string_pass = mod_base_string_pass
	new_qual_string_pass = mod_qual_string_pass
	new_map_string_pass = mod_map_string_pass
	for key in cdict.keys():
		if len(new_base_string + cdict[key][0]) != len(new_qual_string + cdict[key][1]):
			print("error in function:", string, qual_string)
			return []
		else:
			new_base_string += cdict[key][0]
			new_qual_string += cdict[key][1]
			new_map_string += cdict[key][2]
			new_base_string_pass += cdict[key][3]
			new_qual_string_pass += cdict[key][4]
			new_map_string_pass += cdict[key][5]
	return ([new_base_string, new_qual_string, new_map_string, new_base_string_pass, new_qual_string_pass, new_map_string_pass], insert, delete, len(new_base_string) + insert + delete)



def count_variant_bases_in_string(string):
	N = 0
	for nucleotide in "ACTGactg":
		N += string.count(nucleotide)
	return N



def pileup_to_rebase_frequent_variant_whole_file(input, output, indel, depth): #, normal_bedcov_list, tumor_bedcov_list):
# input: a pileup file
# output: a rebase file
# use rebase_with_frequent_varient_allele to convert the whole pileup file

	global GERMLINE_VARIANT_COUNT 
	GERMLINE_VARIANT_COUNT = get_germline_variant_count_threshold(depth)
	global NORMAL_COUNT_BINOM 
	NORMAL_COUNT_BINOM = get_normal_count_binom_threshold(depth)

	f = open(input, 'r')
	g = open(output, 'w')
	indel_file = open(indel, 'w')
#	h = open(output + "_lowcount", 'w')
	current = 0
	for line in f:
		sp = line.strip().split('\t')
		#if sp[0] == 'chr':
		#	continue
		if sp[3] == '0':
			continue
		if sp[4][0] == "\"":
			sp[4] = sp[4][1:(len(sp[4]) -1)]
		rebase_tumor = rebase_with_frequent_variant_allele(sp[4], sp[5], sp[6])
		tumor_insert = rebase_tumor[1]
		tumor_delete = rebase_tumor[2]
		tumor_count = rebase_tumor[3]
		rebase_tumor = rebase_tumor[0]
		if sp[7] != '0':
			rebase_normal = rebase_with_frequent_variant_allele(sp[8], sp[9], sp[10])
			normal_insert = rebase_normal[1]
			normal_delete = rebase_normal[2]
			normal_count = rebase_normal[3]
			rebase_normal = rebase_normal[0]
		else:
			continue
		if sp[11] != '0':
			rebase_extendedFrags = rebase_with_frequent_variant_allele(sp[12],sp[13],sp[14])
			extendedFrags_insert = rebase_extendedFrags[1]
			extendedFrags_delete = rebase_extendedFrags[2]
			extendedFrags_count = rebase_extendedFrags[3]
			rebase_extendedFrags = rebase_extendedFrags[0]
		else:
			rebase_extendedFrags = ["", "", "", "", "", ""]
			extendedFrags_insert = 0
			extendedFrags_delete = 0
			extendedFrags_count = 0
		if sp[15] != '0':
			rebase_notCombined = rebase_with_frequent_variant_allele(sp[16],sp[17],sp[18])
			notCombined_insert = rebase_notCombined[1]
			notCombined_delete = rebase_notCombined[2]
			notCombined_count = rebase_notCombined[3]
			rebase_notCombined = rebase_notCombined[0]
		else:
			rebase_notCombined = ["", "", "", "", "", ""]
			notCombined_insert = 0
			notCombined_delete = 0
			notCombined_count = 0
		if normal_delete >=2 or tumor_delete >= 2 or extendedFrags_delete >= 2 or notCombined_delete >= 2:
			indel_file.write(sp[0] + '\t' + sp[1] + '\t' + "delete" + '\t' + str(normal_delete) + '\t' + str(normal_count) + '\t' + str(tumor_delete) + '\t' + str(tumor_count) + '\t' + str(extendedFrags_delete) + '\t' + str(extendedFrags_count) + '\t' + str(notCombined_delete) + '\t' + str(notCombined_count) + '\n')
			continue
		if normal_insert >= 2 or tumor_insert >= 2 or extendedFrags_insert >= 2 or notCombined_insert >=2:
			indel_file.write(sp[0] + '\t' + sp[1] + '\t' + "insert" + '\t' + str(normal_insert) + '\t' + str(normal_count) + '\t' + str(tumor_insert) + '\t' + str(tumor_count) + '\t' + str(extendedFrags_insert) + '\t' + str(extendedFrags_count) + '\t' + str(notCombined_insert) + '\t' + str(notCombined_count) + '\n')
			continue
		if rebase_tumor == [] or rebase_normal == []:
			continue
		else:
			chr = sp[0]
			pos = int(sp[1])
#			normal_depth, tumor_depth, current = find_coverage_in_exon(chr, pos, normal_bedcov_list, tumor_bedcov_list, current)
			if len(rebase_tumor[0]) != len(rebase_tumor[1]) or len(rebase_normal[0]) != len(rebase_normal[1]) or len(rebase_extendedFrags[0]) != len(rebase_extendedFrags[1]) or len(rebase_notCombined[0]) != len(rebase_notCombined[1]):
				print("error:", i, rebase_tumor, rebase_normal)
				break 
#			basestring, quallist, maplist, qualstring, mapstring = filter_low_qual_with_string(rebase_tumor[0], string_to_qual(rebase_tumor[1]), string_to_qual(rebase_tumor[2]), rebase_tumor[1], rebase_tumor[2])
#			basestring_normal, quallist_normal, maplist_normal, qualstring_normal, mapstring_normal = filter_low_qual_with_string(rebase_normal[0], string_to_qual(rebase_normal[1]), string_to_qual(rebase_normal[2]), rebase_normal[1], rebase_normal[2])
#			basestring_extendedFrags, quallist_extendedFrags, maplist_extendedFrags, qualstring_extendedFrags, mapstring_extendedFrags = filter_low_qual_with_string(rebase_extendedFrags[0], string_to_qual(rebase_extendedFrags[1]), string_to_qual(rebase_extendedFrags[2]), rebase_extendedFrags[1], rebase_extendedFrags[2])
#			basestring_notCombined, quallist_notCombined, maplist_notCombined, qualstring_notCombined, mapstring_notCombined = filter_low_qual_with_string(rebase_notCombined[0], string_to_qual(rebase_notCombined[1]), string_to_qual(rebase_notCombined[2]), rebase_notCombined[1], rebase_notCombined[2])
			basestring = rebase_tumor[3]
			basestring_normal = rebase_normal[3]
			basestring_normal_all = rebase_normal[0]
			basestring_extendedFrags = rebase_extendedFrags[3]
			basestring_notCombined = rebase_notCombined[3]
			basestring_extendedFrags_all = rebase_extendedFrags[0]
			basestring_notCombined_all = rebase_notCombined[0]
			if count_variant_bases_in_string(basestring) == 0 and count_variant_bases_in_string(basestring_extendedFrags) + count_variant_bases_in_string(basestring_notCombined) == 0 :
				continue
			if len(basestring) == 0 or len(basestring_normal) == 0:
				continue
			basecount_notCombined = count_base(basestring_notCombined)
			basecount_extendedFrags = count_base(basestring_extendedFrags)
			basecount_notCombined_all = count_base(basestring_notCombined_all)
			basecount_extendedFrags_all = count_base(basestring_extendedFrags_all)
			variant_base_all = find_major_variant(basecount_notCombined_all, basecount_extendedFrags_all)
			variant_base = find_major_variant(basecount_notCombined, basecount_extendedFrags)
			basecount_tumor_all = count_base(rebase_tumor[0])
			basecount_tumor = count_base(basestring)
			triallele_tumor_all = filter_triallelic_position(basecount_tumor_all, variant_base_all, depth)
			triallele_tumor = filter_triallelic_position(basecount_tumor, variant_base, depth)
			cnt_var_normal = basestring_normal.upper().count(variant_base)
			cnt_alt_normal = sum([1 for i in basestring_normal.upper() if i != 'R'])
			cnt_normal = len(basestring_normal)
			cnt_var_tumor = basestring.upper().count(variant_base)
			basestring_merge = basestring_extendedFrags + basestring_notCombined
			if len(basestring_merge) < DEPTH_FOR_DETECTION_TUMOR or len(basestring) < DEPTH_FOR_DETECTION_TUMOR:
				continue
			if len(basestring_normal) < DEPTH_FOR_DETECTION_NORMAL:
				continue
			VAF_normal = float(cnt_var_normal)/float(len(basestring_normal))
			cnt_var_tumor_merge = basestring_merge.upper().count(variant_base)
			VAF_tumor_merge = float(cnt_var_tumor_merge)/float(len(basestring_merge))
			VAF_tumor = float(cnt_var_tumor)/float(len(basestring))
			if VAF_normal > SOMATIC_VAF_THRESHOLD_IN_NORMAL:
				continue
			elif cnt_var_normal > NORMAL_COUNT_VAR or cnt_alt_normal > NORMAL_COUNT_ALT:
				continue
			elif basestring_normal_all.upper().count(variant_base_all) > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER or  basestring_normal.upper().count(variant_base) > GERMLINE_VARIANT_COUNT or basestring_normal_all.upper().count(variant_base) > GERMLINE_VARIANT_COUNT_BEFORE_QUAL_FILTER:
				continue 
			elif triallele_tumor_all == "F" or triallele_tumor == "F":
				continue
			elif variant_base != variant_base_all:
				continue
			elif VAF_normal > SOMATIC_VAF_THRESHOLD_IN_NORMAL:
				continue
			elif cnt_var_normal > 0 and binom.cdf(cnt_var_normal, len(basestring_normal), VAF_tumor) > NORMAL_COUNT_BINOM:
				continue
			elif cnt_var_normal > NORMAL_COUNT_VAR or cnt_alt_normal > NORMAL_COUNT_ALT:
				continue
			elif cnt_var_normal > 0 and (VAF_normal >= VAF_tumor_merge or VAF_normal >= VAF_tumor or cnt_var_normal >= cnt_var_tumor or cnt_var_normal >= cnt_var_tumor_merge):
				continue
			elif count_variant_bases_in_string(basestring) >= COUNT_VAR_TUMOR_HIGHQUAL or count_variant_bases_in_string(basestring_extendedFrags) + count_variant_bases_in_string(basestring_notCombined) >= COUNT_VAR_TUMOR_HIGHQUAL :
				 g.write('\t'.join(sp[0:4])+'\t'+'\t'.join(rebase_tumor[0:3])+'\t'+sp[7]+'\t'+'\t'.join(rebase_normal[0:3])+'\t'+sp[11]+'\t'+'\t'.join(rebase_extendedFrags[0:3])+'\t'+sp[15]+'\t'+'\t'.join(rebase_notCombined[0:3])+'\t'+'\t'.join(rebase_tumor[3:6])+'\t'+'\t'.join(rebase_normal[3:6])+'\t'+'\t'.join(rebase_extendedFrags[3:6])+'\t'+'\t'.join(rebase_notCombined[3:6])+'\n')
	indel_file.close()
	g.close()
#	h.close()
	f.close()


def find_coverage_in_exon(chr, pos, normal_bedcov_list, tumor_bedcov_list, current):
	index_all = range(current, len(normal_bedcov_list[0]))
	for i in index_all:
		if normal_bedcov_list[0][i] == chr and normal_bedcov_list[1][i] > pos:
			flag = 0
			break
		if normal_bedcov_list[0][i] == chr and normal_bedcov_list[1][i] <= pos and normal_bedcov_list[2][i] >= pos:
			flag = 1
			break
		if normal_bedcov_list[0][i] == chr and normal_bedcov_list[2][i] < pos and i == len(normal_bedcov_list[0]) - 1:
			flag = 0
			break
	current = i
	if flag == 1:
		print(i, current, len(normal_bedcov_list[3]))
		normal_depth = normal_bedcov_list[3][i]
		tumor_depth = tumor_bedcov_list[3][i]
	elif flag == 0:
		normal_depth = "not_for_CN"
		tumor_depth = "not_for_CN"
	return normal_depth, tumor_depth, current

def read_bedcov_to_list(fbedcov):
	bedcov = open(fbedcov, 'r')
	bedcov_chr = []
	bedcov_start = []
	bedcov_end = []
	bedcov_coverage = []
	for line in bedcov:
		sp = line.strip().split('\t')
		bedcov_chr.append(sp[0])
		bedcov_start.append(int(sp[1]))
		bedcov_end.append(int(sp[2]))
		bedcov_coverage.append(float(sp[len(sp) - 1]))
	bedcov.close()
	bedcov_list = [bedcov_chr, bedcov_start, bedcov_end, bedcov_coverage]
	return bedcov_list

#LARGE_NUMBER_FOR_NORMALIZATION
def adjust_bedcov_coverage_with_flagstat(bedcov_list, fflagstat):
	flagstat = open(fflagstat, 'r')
	tmp = flagstat.readline()
	while tmp != '':
		sp = tmp.split(' ')
		if sp[3] == "mapped":
			break
		else:
			tmp = flagstat.readline()
	read_count = float(sp[0])
	flagstat.close()
	#print tmp
	#print read_count
	#bedcov_coverage = [ i * LARGE_NUMBER_FOR_NORMALIZATION / read_count for i in bedcov_list[3] ]
	bedcov_coverage = [ float(bedcov_list[3][i])/float(bedcov_list[2][i] - bedcov_list[1][i]) for i in range(len(bedcov_list[3])) ]
	avg_coverage = sum(bedcov_coverage)/float(len(bedcov_coverage))
	std_coverage = std(bedcov_coverage, dtype=np.float64)
	bedcov_coverage = [ (i-avg_coverage)/std_coverage for i in bedcov_coverage]
	bedcov_list[3] = bedcov_coverage
	return

def extract_coverage_in_exon(fnormal_bedcov, fnormal_flagstat, ftumor_bedcov, ftumor_flagstat):
	normal_bedcov_list = read_bedcov_to_list(fnormal_bedcov)
	adjust_bedcov_coverage_with_flagstat(normal_bedcov_list, fnormal_flagstat)
	tumor_bedcov_list = read_bedcov_to_list(ftumor_bedcov)
	adjust_bedcov_coverage_with_flagstat(tumor_bedcov_list, ftumor_flagstat)
	return normal_bedcov_list, tumor_bedcov_list

if __name__ == "__main__":
	finput = sys.argv[1]
	foutput = sys.argv[2]
	findel = sys.argv[3]
	depth = float(sys.argv[4])

#	fnormal_flagstat = sys.argv[3]
#	fnormal_bedcov = sys.argv[4]
#	ftumor_flagstat = sys.argv[5]
#	ftumor_bedcov = sys.argv[6]
#	normal_bedcov_list, tumor_bedcov_list = extract_coverage_in_exon(fnormal_bedcov, fnormal_flagstat, ftumor_bedcov, ftumor_flagstat)
	pileup_to_rebase_frequent_variant_whole_file(finput, foutput, findel, depth)#, normal_bedcov_list, tumor_bedcov_list)


