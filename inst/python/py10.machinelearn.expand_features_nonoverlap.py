import sys

WINDOW_SIZE = 3
HOMOPOLYMER_SIZE = 5
CIGAR_LETTERS = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
REF_LETTERS = ['A', 'C', 'G', 'T']
VAR_LETTERS = ['A', 'C', 'G', 'T', 'N', 'D']
file_id_list = sys.argv[1].split(',')
file_list = sys.argv[2].split(',')

# tab-delimited text file format
#
# an example:
# 0       100M    22      29      29      28      28      29      29      G|G     T|T     C|C     G|A     C|C     C|C     C|CMM       M       M       M       M       M       -1      0       281     1       0       1       163     83      10      100MNA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA NA       NA      NA      NA      NA      -1      NA      chr1    1461917 chr1    1462098 chr1    1461963 GTCGCCC
# 27      100M    25      28      30      30      29      31      30      G|G     T|T     C|C     G|A     C|C     C|C     C|CMM       M       M       M       M       M       -1      0       168     1       0       1       99      147     27      100MNA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA NA       NA      NA      NA      NA      -1      NA      chr1    1461900 chr1    1461968 chr1    1461963 GTCGCCC
# 33      100M    25      30      31      32      31      31      32      G|G     T|T     C|C     G|A     C|C     C|C     C|CMM       M       M       M       M       M       -1      0       185     1       0       1       163     83      33      100MNA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA NA       NA      NA      NA      NA      -1      NA      chr1    1461887 chr1    1461972 chr1    1461963 GTCGCCC
#
# mapping quality R1
# cigar string R1
# base quality R1 -3
# base quality R1 -2
# base quality R1 -1
# base quality R1 0
# base quality R1 1
# base quality R1 2
# base quality R1 3
# ref:read R1 -3   reference base and base on the read separated by a colon (:) at position -3 (the third base left from the possible variant)
# ref:read R1 -2
# ref:read R1 -1
# ref:read R1 0
# ref:read R1 1
# ref:read R1 2
# ref:read R1 3
# distance to nearby indel on read R1: -1 no indel
# if it is in homopolymer greater than 5 R1
#
# fragment size
# if both reads are on the same chromosome
# strand direction R1
# starnd direction R2
# flag R1
# flag R2
#
# mapping quality R2
# cigar string R2
# base quality R2 -3
# base quality R2 -2
# base quality R2 -1
# base quality R2 0
# base quality R2 1
# base quality R2 2
# base quality R2 3
# ref:read R2 -3   reference base and base on the read separated by a colon (:) at position -3 (the third base left from the possible variant)
# ref:read R2 -2
# ref:read R2 -1
# ref:read R2 0
# ref:read R2 1
# ref:read R2 2
# ref:read R2 3
# distance to nearby indel on read R2: -1 no indel
# if it is in homopolymer greater than 5 R2
#
# chromosome R1
# mapping position on reference R1
# chromosome R2
# mapping position on reference R2
# variant chromosome 
# variant position on reference 
# sequencing context around variant window = 7


def transform_CIGAR_to_expanded_features(unexpanded_features): 
	# transform data line by line
	# unexpanded_features is a CIGAR string to be expanded
	# unexpanded_features is finally changed to a set of categorical data, where each column is 1 or 0 indicating whether a character is in the CIGAR string or not
	expanded_features = []
	for i in CIGAR_LETTERS:
		if i in unexpanded_features:
			expanded_features.append('1')
		else:
			expanded_features.append('0')
	return expanded_features


def expand_seq_feature(unexpand_string):
	expanded_ref = ['0'] * (len(REF_LETTERS))
	expanded_var = ['0'] * (len(VAR_LETTERS))
	if unexpand_string == 'NA':
		return expanded_ref	+ expanded_var
	ref, var = unexpand_string.split('|')
#	print ref, var
	expanded_ref[REF_LETTERS.index(ref)] = '1'
	expanded_var[VAR_LETTERS.index(var)] = '1'
	return expanded_ref	+ expanded_var


CIGAR_CODES = ['D', 'M', 'N', 'X', '=']

def expand_CIGAR_code(cigar_code):
        expanded_cigar_code = ['0'] * (len(CIGAR_CODES))
        if cigar_code == 'NA':
                return expanded_cigar_code
        expanded_cigar_code[CIGAR_CODES.index(cigar_code)] = '1'
        return expanded_cigar_code


def transform_seq_features_to_expanded_features(unexpanded_features):
	expanded_features = []
	for unexpanded_string in unexpanded_features:
		expanded_string = expand_seq_feature(unexpanded_string.upper())
		expanded_features += expanded_string
	return expanded_features


def transform_CIGAR_codes_to_expanded_features(cigar_code_list):
        expanded_features = []
        for cigar_code in cigar_code_list:
                expanded_string = expand_CIGAR_code(cigar_code)
                expanded_features += expanded_string
        return expanded_features


def generate_header():
	header_info = ["mapping_qual"]
	header_info += ["cigar:" + i for i in CIGAR_LETTERS]
	header_info += ["base_qual_" + str(i - WINDOW_SIZE) for i in range(2*WINDOW_SIZE + 1)]
	for i in range(2*WINDOW_SIZE + 1):
		header_info += ["ref_" + str(i - WINDOW_SIZE) + ":" + j for j in REF_LETTERS]
		header_info += ["var_" + str(i - WINDOW_SIZE) + ":" + j for j in VAR_LETTERS]
	for i in range(2*WINDOW_SIZE + 1):
		header_info += ["cigar_" + str(i - WINDOW_SIZE) + ":" + j for j in CIGAR_CODES ]
	header_info += ["dist_indel", "homopolymer_read_>"+str(HOMOPOLYMER_SIZE), "strand", "flag"]
	header_relation = ["frag_size", "if_same_chrom", "homopolymer_reference_>"+str(HOMOPOLYMER_SIZE)]
#	header = [i + "_R1" for i in header_info] + header_relation + [i + "_R2" for i in header_info]
	header = [i + "_R1" for i in header_info] + header_relation + ["mapping_qual_R2"] + ["cigar:" + i + "_R2" for i in CIGAR_LETTERS] + ["dist_indel_R2", "homopolymer_>5", "strand_R2", "flag_R2"]
	return header

def transform_line_to_expanded_features(line):
        sp = line.strip().split('\t')
        expanded_features = [sp[0]] #
        expanded_features += transform_CIGAR_to_expanded_features(sp[1]) #
        expanded_features += ['0' if i == "NA" or i == "D"  else i for i in sp[2:(2*WINDOW_SIZE + 3)]] #
        expanded_features += transform_seq_features_to_expanded_features(sp[(2*WINDOW_SIZE + 3):(4*WINDOW_SIZE+4)]) #
        expanded_features += transform_CIGAR_codes_to_expanded_features(sp[(4*WINDOW_SIZE+4):(6*WINDOW_SIZE+5)]) #
        for i in [6*WINDOW_SIZE+5, 6*WINDOW_SIZE + 6, 6*WINDOW_SIZE + 7, 6*WINDOW_SIZE + 8, 6*WINDOW_SIZE + 9,6*WINDOW_SIZE + 10, 6*WINDOW_SIZE + 11, 6*WINDOW_SIZE + 12, 12*WINDOW_SIZE+18, 12*WINDOW_SIZE+19, 12*WINDOW_SIZE + 20]: # [16,17,18,19,20,21,22,23,40,41]:
                if sp[i] == "NA":
                        sp[i] = "-1"
        expanded_features += [sp[6*WINDOW_SIZE+5], sp[6*WINDOW_SIZE+6], sp[6*WINDOW_SIZE+9], sp[6*WINDOW_SIZE+11], sp[6*WINDOW_SIZE+7], sp[6*WINDOW_SIZE+8], sp[12*WINDOW_SIZE + 20], sp[6*WINDOW_SIZE + 13]] # [sp[16], sp[17], sp[20], sp[22], sp[18], sp[19], sp[24]]
        expanded_features += transform_CIGAR_to_expanded_features(sp[6*WINDOW_SIZE + 14]) # (sp[25])
#        expanded_features += ['0' if i == "NA" or i == "D" else i for i in sp[33:40]] # sp[26:33]]
#        expanded_features += transform_seq_features_to_expanded_features(sp[40:47]) #(sp[33:40])
#        expanded_features += transform_CIGAR_codes_to_expanded_features(sp[47:54])
        expanded_features += [sp[12*WINDOW_SIZE + 18], sp[12*WINDOW_SIZE + 19], sp[6*WINDOW_SIZE+10], sp[6*WINDOW_SIZE+12]] #[sp[40], sp[41], sp[21], sp[23]]
        return expanded_features
#	sp = line.strip().split('\t')
#	expanded_features = [sp[0]]
#	expanded_features += transform_CIGAR_to_expanded_features(sp[1])
#	expanded_features += ['0' if i == "NA" or i == "D" else i for i in sp[2:9]]
#	expanded_features += transform_seq_features_to_expanded_features(sp[9:16])
#	for i in [16,17,18,19,29,21,22,23,40,41]:
#		if sp[i] == "NA":
#			sp[i] = "-1"
#	expanded_features += [sp[16], sp[17], sp[20], sp[22], sp[18], sp[19], sp[24]]
#	expanded_features += transform_CIGAR_to_expanded_features(sp[25])
##	expanded_features += ['0' if i == "NA" else i for i in sp[26:33]]
##	expanded_features += transform_seq_features_to_expanded_features(sp[33:40])
#	expanded_features += [sp[40], sp[41], sp[21], sp[23]]	
#	return expanded_features


for i in range(len(file_list)):
	file = file_list[i]
	ind = file_id_list[i]
	with open(file) as f:
		f_expand_list = []
		f_index_list = []
		n = 0
		header = generate_header()
		for line in f:
			n += 1 # record index is 1-based
			f_index_list.append(ind + '_' + str(n))
			f_expand_list.append(transform_line_to_expanded_features(line))
	with open(file + ".expand", 'w') as g:
		g.write("id" + '\t' + '\t'.join(header) + '\n')
		for j in range(n):
			g.write(f_index_list[j] + '\t' + '\t'.join(f_expand_list[j]) + '\n')



