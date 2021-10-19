from parameter import *
from _filter import *
from _jenks import *
from _probability import *
import numpy as np
import sys


def generate_intermediate_result_one_line(line):
	in_VAF = 'F'
	sp = line.strip().split('\t')
	sp = [i[1:(len(i)-1)] for i in sp]
	chrom = sp[0]
	pos = sp[1]
	#print chrom, pos
	variant_base = sp[9]
	VAF_merge = float(sp[14])
	VAF = float(sp[11])
	joint_genotype = sp[10]
	joint_genotype_merge = sp[13]
	basestring = sp[16]
	quallist = string_to_qual(sp[17])
	maplist = string_to_qual(sp[18])
	basestring_extendedFrags = sp[22]
	quallist_extendedFrags = string_to_qual(sp[23])
	maplist_extendedFrags = string_to_qual(sp[24])
	basestring_notCombined = sp[25]
	quallist_notCombined = string_to_qual(sp[26])
	maplist_notCombined = string_to_qual(sp[27])
	cnt_normal = len(sp[19])
	basestring_merge = basestring_extendedFrags + basestring_notCombined
	quallist_merge = np.concatenate([quallist_extendedFrags, quallist_notCombined])
	maplist_merge = np.concatenate([maplist_extendedFrags, maplist_notCombined])
	cnt_var_tumor = int(sp[6])
	cnt_var_normal = int(sp[4])
	cnt_var_tumor_merge = int(sp[8])
	basecount_extendedFrags = count_base(basestring_extendedFrags)
	basecount_notCombined = count_base(basestring_notCombined)
	basecount = count_base(basestring)
	strand_bias_merge = filter_strand_bias_merge(basecount_notCombined, basecount_extendedFrags,variant_base)
	both_strand_merge, above_average_merge = filter_both_strand_above_average_merge(basecount_notCombined, basecount_extendedFrags)
	#print basestring, basecount
	strand_bias_unmerge = filter_strand_bias_unmerge(basecount, variant_base)
	both_strand_unmerge, above_average_unmerge = filter_both_strand_above_average_unmerge(basecount)
	all_mapping_quality_unmerge = filter_all_mapping_quality(maplist)
	all_mapping_quality_merge = filter_all_mapping_quality(maplist_merge)
	variant_mapping_quality_unmerge = filter_variant_mapping_quality(maplist, basestring, variant_base)
	variant_mapping_quality_merge = filter_variant_mapping_quality(maplist_merge, basestring_merge, variant_base)
	variant_base_quality_unmerge = filter_variant_base_quality(quallist, basestring, variant_base)
	variant_base_quality_merge = filter_variant_base_quality(quallist_merge, basestring_merge, variant_base)
	not_in_dbSNP = filter_dbSNP(chrom, pos, variant_base)
	supporting_read_unmerge = filter_supporting_reads_unmerge(basestring, variant_base)
	supporting_frag_merge = filter_supporting_fragment_merge(basestring_merge, variant_base)
	VAF_filter_unmerge, VAF_descript_unmerge = filter_VAF(VAF, depth)
	VAF_filter_merge, VAF_descript_merge = filter_VAF(VAF_merge, depth)
	cnt = len(basestring)
	cnt_merge = len(basestring_merge)
	normal_coverage = filter_normal_coverage(cnt_normal)
	tumor_coverage_unmerge, coverage_descript_unmerge = filter_tumor_coverage(cnt)
	tumor_coverage_merge, coverage_descript_merge = filter_tumor_coverage(cnt_merge)
	if normal_coverage == 'T':
		if [strand_bias_merge, both_strand_merge, above_average_merge, all_mapping_quality_merge, variant_mapping_quality_merge, variant_base_quality_merge, supporting_frag_merge, VAF_filter_merge] == ['T'] * 8 and cnt_merge > 50:
			if [strand_bias_unmerge, both_strand_unmerge, above_average_unmerge, all_mapping_quality_unmerge, variant_mapping_quality_unmerge, variant_base_quality_unmerge, supporting_read_unmerge, VAF_filter_unmerge] == ['T'] * 8 and cnt > 50 :
				in_VAF = 'T'
	output = sp[0:3] + [sp[9]] + sp[3:7] + sp[10:13] + sp[7:9] + sp[13:16] + sp[16:28] + [normal_coverage]
	output += [strand_bias_merge, both_strand_merge, above_average_merge, all_mapping_quality_merge, variant_mapping_quality_merge, variant_base_quality_merge, supporting_frag_merge, VAF_filter_merge, VAF_descript_merge, tumor_coverage_merge, coverage_descript_merge]
	output += [strand_bias_unmerge, both_strand_unmerge, above_average_unmerge, all_mapping_quality_unmerge, variant_mapping_quality_unmerge, variant_base_quality_unmerge, supporting_read_unmerge, VAF_filter_unmerge, VAF_descript_unmerge, tumor_coverage_unmerge, coverage_descript_unmerge]
	output = list(map(str, output))
	return output, in_VAF, VAF_merge	


def generate_intermediate_result_whole_file(filename, output, database):
# input: name of the variant calling result, name of the filtration result
# output: none
	import_dbSNP(database)
	f = open(filename, 'r')
	g = open(output, 'w')
	VAF_list = []
	for line in f:
		output, in_VAF, VAF_merge = generate_intermediate_result_one_line(line)
		if in_VAF == 'T':
			VAF_list.append(VAF_merge)
		g.write('\t'.join(output) + '\n')
	g.close()
	f.close()
	return VAF_list


def header_intermediate_file():
	header_list = ["chromosome","position","reference","variant", "cnt_normal", "cnt_var_normal", "cnt_tumor", "cnt_var_tumor","joint_genotype_MAP","VAF_tumor","confidence", "cnt_tumor_merge", "cnt_var_tumor_merge","joint_genotype_MAP_merge","VAF_tumor_merge","confidence_merge","basestring","qual_string","map_string","basestring_normal","qual_string_normal","map_string_normal","basestring_ef","qual_string_ef","map_string_ef","basestring_nc","qual_string_nc","mapstring_nc","normal_coverage","strand_bias_merge","both_strand_merge","above_average_merge","all_mapping_quality_merge","variant_mapping_quality_merge","variant_base_quality_merge","supporting_frag_merge","VAF_filter_merge","VAF_descript_merge","tumor_coverage_merge","descript_merge","strand_bias_unmerge","both_strand_unmerge","above_average_unmerge","all_mapping_quality_unmerge","variant_mapping_quality_unmerge","variant_base_quality_unmerge","supporting_frag_unmerge","VAF_filter_unmerge","VAF_descript_unmerge","tumor_coverage_unmerge","descript_unmerge"]
	header = '\t'.join(header_list)
	return header


def header_record_file():
	header_list = ["chromosome","position","reference","variant", "cnt_normal", "cnt_var_normal", "cnt_tumor", "cnt_var_tumor","joint_genotype_MAP","VAF_tumor","confidence","cnt_tumor_merge","cnt_var_tumor_merge","joint_genotype_MAP_merge","VAF_tumor_merge","confidence_merge","basestring","qual_string","map_string","basestring_normal","qual_string_normal","map_string_normal","basestring_ef","qual_string_ef","map_string_ef","basestring_nc","qual_string_nc","mapstring_nc","normal_coverage","strand_bias_merge","both_strand_merge","above_average_merge","all_mapping_quality_merge","variant_mapping_quality_merge","variant_base_quality_merge","supporting_frag_merge","VAF_filter_merge","VAF_descript_merge","tumor_coverage_merge","descript_merge","binomial_test_merge","status_merge","strand_bias_unmerge","both_strand_unmerge","above_average_unmerge","all_mapping_quality_unmerge","variant_mapping_quality_unmerge","variant_base_quality_unmerge","supporting_frag_unmerge","VAF_filter_unmerge","VAF_descript_unmerge","tumor_coverage_unmerge","descript_unmerge","binomial_test_unmerge","status_unmerge","decision"]
	header = '\t'.join(header_list)
	return header


def header_output_file():
	header_list = ["chromosome","position","reference","variant","cnt_normal","cnt_var_normal","cnt_tumor","cnt_var_tumor","joint_genotype_MAP","VAF_tumor","confidence","cnt_tumor_merge", "cnt_var_tumor_merge","joint_genotype_MAP_merge","VAF_tumor_merge","confidence_merge"]
	header = '\t'.join(header_list)
	return header



def generate_record_one_line(line):
	sp = line.strip().split('\t')
	chrom = sp[0]
	pos = sp[1]
	ref = sp[2]
	variant_base = sp[3]
	basestring = sp[16]
	basestring_extendedFrags = sp[22]
	basestring_notCombined = sp[25]
	basestring_merge = basestring_extendedFrags + basestring_notCombined
	binom_test_unmerge = filter_supporting_count_by_binomial_test(jenks_estimate, basestring, variant_base)
	binom_test_merge = filter_supporting_count_by_binomial_test(jenks_estimate, basestring_merge, variant_base)
	sp.insert(40, binom_test_merge)
	sp.append(binom_test_unmerge)
	sp = list(map(str, sp))
	return sp



def decide_output_status_one_line(record):
	merge = record[29:37] + [record[38], record[40]]
	#print merge
	unmerge = record[41:49] + [record[50], record[52]]
	#print unmerge
	output = record[0:16]
	if record[28] == "T" and record[41] == "T" and merge[0] == "T" and merge[2:10] == ["T"]*8:
		merge_status = "pass"
	elif record[28] == "T" and record[41] == "T" and merge[0] == "T" and merge[2:5] == ["T"]*4 and merge[9] == "T":
		merge_status = "hold"
	else:
		merge_status = "reject"
	if record[28] == "T" and unmerge[0] == "T" and unmerge[2:10] == ["T"]*8:
		unmerge_status = "pass"
	elif record[28] == "T" and unmerge[0] == "T" and unmerge[2:5] == ["T"]*4 and unmerge[9] == "T":
		unmerge_status = "hold"
	else:
		unmerge_status = "reject"
	#print merge_status, unmerge_status
	if merge_status == "pass" and unmerge_status == "pass":
		return output, merge_status, unmerge_status, "pass"
	elif merge_status == "pass" and unmerge_status == "hold":
		return output, merge_status, unmerge_status, "hold"
	elif merge_status == "hold" and unmerge_status == "pass":
		return output, merge_status, unmerge_status, "hold"
	elif merge_status == "hold" and unmerge_status == "hold":
		return output, merge_status, unmerge_status, "hold"
	elif merge_status == "pass" and unmerge_status == "reject":
		return output, merge_status, unmerge_status, "hold"
	elif merge_status == "reject" and unmerge_status == "pass":
		return output, merge_status, unmerge_status, "hold"
	elif merge_status == "reject" and unmerge_status == "reject":
		return output, merge_status, unmerge_status, "reject"
	elif merge_status == "hold" and unmerge_statis == "reject":
		return output, merge_status, unmerge_status, "hold"
	elif merge_status == "reject" and unmerge_status == "hold":
		return output, merge_status, unmerge_status, "hold"

def generate_record_and_output_whole_file(fintermediate, foutput_pass, foutput_check, frecord, jenks_estimate):
	f = open(fintermediate, 'r')
	r = open(frecord, 'w')
	g = open(foutput_pass, 'w')
	h = open(foutput_check, 'w')
	r_header = header_record_file()
	r.write(r_header + '\n')
	n = 0
	for line in f:
		record = generate_record_one_line(line) ##
		output, merge_status, unmerge_status, status = decide_output_status_one_line(record)
		if status == "pass":
			n += 1
			g.write('\t'.join(output) + '\n')
		elif status == "hold":
			n += 1
			h.write('\t'.join(output) + '\n')	
		record.insert(41, merge_status)
		record.append(unmerge_status)
		record.append(status)
		record_out = ["\"" + i + "\"" for i in record]
		r.write('\t'.join(record_out) + '\n')
	return n


if __name__ == "__main__":
	finput = sys.argv[1]
	fintermediate = sys.argv[2]
	foutput_pass = sys.argv[3]
	foutput_check = sys.argv[4]
	frecord = sys.argv[5]
	fTF = sys.argv[6]
	MERGED_VAF_THRESHOLD = sys.argv[7]
	SNP_database = sys.argv[8]
	depth = float(sys.argv[9])
	VAF_list = generate_intermediate_result_whole_file(finput, fintermediate, SNP_database)
	jenks_estimate, include_number, include_group = final_estimation_with_Jenks(VAF_list)
	n = generate_record_and_output_whole_file(fintermediate, foutput_pass, foutput_check, frecord, jenks_estimate)
	with open(fTF, 'w') as TF:
		TF.write(str(n) + '\t' + str(jenks_estimate) + '\t' + str(format(jenks_estimate/2.0, '.4f')) + '\t' + str(include_number) + '\t' + str(include_group) + '\t' + MERGED_VAF_THRESHOLD + '\n')


