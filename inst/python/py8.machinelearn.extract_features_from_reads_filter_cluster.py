import numpy as np
import sys

loc = []

f = open(sys.argv[1], 'r')

for line in f:
	sp = line.split('\t')
	loc.append( (sp[0], int(sp[2])) )

f.close()



len(loc)

from operator import itemgetter
loc.sort(key = itemgetter(1))
loc.sort(key = itemgetter(0))

n = 0
neighbor_list = []
neighbor = []
flag = 0
for i in range(len(loc)-1):
	if loc[i][0] == loc[i+1][0]:
		if abs(loc[i][1] - loc[i+1][1]) < 20:
			n += 1
			if flag == 0:
				flag = 1
				neighbor += [loc[i], loc[i+1]]
			else:
				neighbor.append(loc[i+1])
#			print loc[i], loc[i+1]
		else:
			if flag == 1:
				neighbor_list.append(neighbor)
			flag = 0
			neighbor = []
	else:
		if flag == 1:
			neighbor_list.append(neighbor)
		flag = 0
		neighbor = []



# In[27]:


from collections import defaultdict

mapping_to_var = defaultdict(lambda: -1)
matrix = []

read_pos_loc = defaultdict(list)

# In[28]:


for i in range(len(neighbor_list)):
    neighbor = neighbor_list[i]
    neighbor_string = [n[0]+"-"+str(n[1]) for n in neighbor]
    matrix.append((neighbor_string, np.zeros((len(neighbor), len(neighbor)))))
    for j in neighbor_string:
        mapping_to_var[j] = i


# In[29]:


WINDOW = 3
HOMOPOLYMER_SIZE = 5
NUCLEOTIDE = ['A', 'C', 'T', 'G']
HOMOPOLYMER = [nucleotide * HOMOPOLYMER_SIZE for nucleotide in NUCLEOTIDE]
BOTH_ADD_CIGAR = "M=X"
REFERENCE_ADD_CIGAR = "DN"
READ_ADD_CIGAR = "IS"

def extract_mapping_quality(read):
	phred = read.split('\t')[4]
	if phred == 255:
		return "NA"
	else:
		return phred


def find_location_on_read(mapping_position, query_position, CIGAR_string, max_length):
	read_current_position = 0
	reference_current_position = mapping_position
	cigar_split = find_all_occurrence_in_string(CIGAR_string)
	cigar_split.insert(0,-1)
	for i in range(0, len(cigar_split) - 1):
		num_bases = int(CIGAR_string[ (cigar_split[i] + 1):cigar_split[i+1] ])
		cigar_case = CIGAR_string[cigar_split[i+1]]
		if cigar_case in BOTH_ADD_CIGAR:
			read_current_position += num_bases
			reference_current_position += num_bases
			if read_current_position > max_length-1 and query_position > reference_current_position:
				#print "exceed max length", read_current_position, max_length, CIGAR_string, query_position, reference_current_position
				return "NA"
			if reference_current_position >= query_position:
				if read_current_position < reference_current_position - query_position:
					#print "exceed left end", reference_current_position, query_position, CIGAR_string
					return "NA"
				else:
					position_on_read = read_current_position - (reference_current_position - query_position)
					if position_on_read > max_length-1 :
						return "NA"
					else:
						return read_current_position - (reference_current_position - query_position)
		elif cigar_case in REFERENCE_ADD_CIGAR:
			reference_current_position += num_bases
			if reference_current_position >= query_position:
				return "D"
		elif cigar_case in READ_ADD_CIGAR:
			read_current_position += num_bases
			if read_current_position > max_length-1 :
				#print "exceed max length", read_current_position, max_length, CIGAR_string
				return "NA"
	#print "ohoh", mapping_position, query_position, CIGAR_string, max_length
	return "NA"


def extract_cigar_base_given_position(mapping_position, query_position, CIGAR_string, max_length):
	read_current_position = 0
	reference_current_position = mapping_position
	cigar_split = find_all_occurrence_in_string(CIGAR_string)
	cigar_split.insert(0,-1)
	for i in range(0, len(cigar_split) - 1):
		num_bases = int(CIGAR_string[ (cigar_split[i] + 1):cigar_split[i+1] ])
		cigar_case = CIGAR_string[cigar_split[i+1]]
		if cigar_case in BOTH_ADD_CIGAR:
			read_current_position += num_bases
			reference_current_position += num_bases
			if read_current_position > max_length-1 and query_position > reference_current_position:
				#print "exceed max length", read_current_position, max_length, CIGAR_string, query_position, reference_current_position
				return "NA"
			if reference_current_position >= query_position:
				if read_current_position < reference_current_position - query_position:
					#print "exceed left end", reference_current_position, query_position, CIGAR_string
					return "NA"
				else:
					position_on_read = read_current_position - (reference_current_position - query_position)
					if position_on_read > max_length-1 :
						return "NA"
					else:
						return cigar_case
		elif cigar_case in REFERENCE_ADD_CIGAR:
			reference_current_position += num_bases
			if reference_current_position >= query_position:
				return cigar_case
		elif cigar_case in READ_ADD_CIGAR:
			read_current_position += num_bases
			if read_current_position > max_length-1 :
				#print "exceed max length", read_current_position, max_length, CIGAR_string
				return "NA"
	#print "ohoh", mapping_position, query_position, CIGAR_string, max_length
	return "NA"



def extract_cigar_base(mapping_position, query_position, CIGAR_string, max_length):
	base = []
	for distance in range(-WINDOW, WINDOW + 1):
		base.append( extract_cigar_base_given_position(mapping_position, query_position + distance, CIGAR_string, max_length) )
	return base


def extract_base_quality_given_distance(variant_position, mapping_position, distance, quality_string, CIGAR_string):
	location_on_read = find_location_on_read(mapping_position, variant_position + distance, CIGAR_string, len(quality_string))
	#print location_on_read, len(quality_string)
	if location_on_read == "NA":
		return "NA"
	elif location_on_read == "D":
		return "D"
	else:
		quality_base = quality_string[location_on_read]
		return str(ord(quality_base)-33)

def extract_base_quality(variant_position, mapping_position, quality_string, CIGAR_string):
	quality = []
	for distance in range(-WINDOW, WINDOW + 1):
		quality.append( extract_base_quality_given_distance(variant_position, mapping_position, distance, quality_string, CIGAR_string) )
	return quality

def extract_reference_base_read_base_given_distance(variant_position, mapping_position, distance, base_string, reference_string, CIGAR_string):
	#print reference_string, distance, WINDOW
	ref = reference_string[distance + WINDOW]
	location_on_read = find_location_on_read(mapping_position, variant_position + distance, CIGAR_string, len(base_string))
	if location_on_read == "NA":
		return "NA"
	elif location_on_read == "D":
		return ref + "|D"
	else:
		onread = base_string[location_on_read]
		return ref + '|' + onread

def extract_reference_base_read_base(variant_position, mapping_position, base_string, reference_string, CIGAR_string):
	bases = []
	for distance in range(-WINDOW, WINDOW + 1):
		bases.append( extract_reference_base_read_base_given_distance(variant_position, mapping_position, distance, base_string, reference_string, CIGAR_string) )
	return bases

CIGAR_CHARACTER = "MSIDHNPX="

def find_all_occurrence_in_string(string):
	return [i for (i, c) in enumerate(string) if c in CIGAR_CHARACTER]

def find_indel_position(cigar_split, cigar_string, position):
	insertion_start = []
	current = position
	string_index = [i for i in cigar_split]
	string_index.insert(0,-1)
	#print cigar_split, string_index
	for i in range(len(cigar_split)):
		#print insertion_start
		id = cigar_split[i]
		#print i, cigar_split[i], string_index[i],string_index[i+1]
		length = int(cigar_string[(string_index[i]+1):string_index[i+1]])
		#print id, length, i
		if cigar_string[id] == "I":
			insertion_start.append(current - 1)
			insertion_start.append(current)
		elif cigar_string[id] == "D":
			start = current - 1
			current += length
			end = current + 1
			insertion_start = insertion_start + list(range(start, end))
		elif cigar_string[id] in "MX=N":
			current += length
	return insertion_start

#???
def extract_nearby_indel(variant_position, mapping_position, CIGAR_string):
	cigar_split = find_all_occurrence_in_string(CIGAR_string)
	all_location_involved_in_indel = find_indel_position(cigar_split, CIGAR_string, mapping_position)
	if all_location_involved_in_indel == []:
		return "-1"
	else:
		distance_to_variant_position = [abs(i - variant_position) for i in all_location_involved_in_indel]
		return str(min(distance_to_variant_position))


def extract_homopolymer_on_read(variant_position, mapping_position, base_string, CIGAR_string):
	location_on_read = find_location_on_read(mapping_position, variant_position, CIGAR_string, len(base_string))
	if location_on_read == "D" or location_on_read == "NA":
		return "NA"
	neighborhood = ""
	for i in range( - HOMOPOLYMER_SIZE + 1, HOMOPOLYMER_SIZE):
		location_on_read = find_location_on_read(mapping_position, variant_position + i, CIGAR_string, len(base_string))
		if location_on_read == "NA" or location_on_read == "D":
			neighborhood += '0'
		else:
			neighborhood += base_string[location_on_read]
	for i in range(HOMOPOLYMER_SIZE):
		if neighborhood[i:(i + HOMOPOLYMER_SIZE)] in HOMOPOLYMER:
			return '1'
	return '0'


def extract_homopolymer_on_reference(ref_string):
	start = max(0, int(0.5 * (len(ref_string) + 1) - HOMOPOLYMER_SIZE))
	end = min(len(ref_string) - HOMOPOLYMER_SIZE + 1, int(0.5 * (len(ref_string) + 1)))
	for i in range(start, end):
		if ref_string[i:(i + HOMOPOLYMER_SIZE)] in HOMOPOLYMER:
			return '1'
	return '0'


def extract_CIGAR_string(read):
	return read.split('\t')[5]



# In[30]:


import sys
from collections import defaultdict
import numpy as np

variant_position_dict = defaultdict(list)
variant_base_dict = defaultdict(list)
context_dict = defaultdict(list)
#with open("SRR6708920.mutect.keep.select_true_variant.bed", 'r') as BED:
with open(sys.argv[1], 'r') as BED:
	for line in BED:
		sp = line.strip().split('\t')
		chr = sp[0]
		pos = sp[2]
		ctx = sp[3]
		ref = sp[4]
		var = sp[5]
		variant_position_dict[chr].append(int(pos))
		variant_base_dict[chr].append(var)
		context_dict[chr].append(ctx.replace('x', ref))

OUTPUT = open(sys.argv[3], 'w')
#with open("SRR6708920.2.select_region_true_variant.variant.paired-reads.recal.sam", 'r') as SAM:
with open(sys.argv[2], 'r') as SAM:
	tmp0 = SAM.readline()
	sp0 = tmp0.strip().split('\t')
	while tmp0 != "":
		tmp1 = tmp0
		sp1 = sp0
		tmp0 = SAM.readline()
		if tmp0 == "":
			break
		sp0 = tmp0.strip().split('\t')
		#print sp0
		if sp0[0] != sp1[0]:
			continue
		chr1 = sp1[2]
		chr0 = sp0[2]
		pos1 = sp1[3]
		pos0 = sp0[3]
		basestring1 = sp1[9]
		basestring0 = sp0[9]
		qualstring1 = sp1[10]
		qualstring0 = sp0[10]
		CIGAR1 = sp1[5]
		CIGAR0 = sp0[5]
		if CIGAR1 == "*" or CIGAR0 == "*":
			continue
		flag1 = sp1[1]
		flag0 = sp0[1]
		binflag1 = bin(int(flag1))
		binflag0 = bin(int(flag0))
		strand1 = binflag1[len(binflag1)-5] # == 1 if its on reverse strand
		strand0 = binflag0[len(binflag0)-5]
		MQ1 = sp1[4]
		MQ0 = sp0[4]
		INFO_1 = [MQ1, CIGAR1]
		INFO_0 = [MQ0, CIGAR0]
		SUP_INFO_1 = [chr1, pos1]
		SUP_INFO_0 = [chr0, pos0]
		if chr1 == chr0:
			size = abs(int(pos0) - int(pos1)) + len(basestring1)
			samechr = 1
			RELATION = [str(size), str(samechr)]
			related_variant = []
			related_read = []
			related_position = []
			for i_var in range(len(variant_position_dict[chr1])):
				var_pos = variant_position_dict[chr1][i_var]
				var_loc_on_read1 = find_location_on_read(int(pos1), int(var_pos), CIGAR1, len(basestring1))
				if var_loc_on_read1 != "NA":
					if var_loc_on_read1 == "D":
#						related_variant.append(i_var)
#						related_read.append(1)
#						related_position.append(var_loc_on_read1)
						continue
					else:
						#print variant_base_dict[chr1][i_var], basestring1[var_loc_on_read1], var_loc_on_read1, basestring1
						if basestring1[var_loc_on_read1] == variant_base_dict[chr1][i_var]:
							related_variant.append(i_var)
							related_read.append(1)
							related_position.append(var_loc_on_read1)
							continue
				var_loc_on_read0 = find_location_on_read(int(pos0), int(var_pos), CIGAR0, len(basestring0))
				if var_loc_on_read0 != "NA":
					if var_loc_on_read0 == "D":
#						related_variant.append(i_var)
#						related_read.append(0)
#						related_position.append(var_loc_on_read0)
						continue
					elif basestring0[var_loc_on_read0] == variant_base_dict[chr0][i_var]:
						related_variant.append(i_var)
						related_read.append(0)
						related_position.append(var_loc_on_read0)
						read_pos_loc[var_pos].append(var_loc_on_read0)
						continue
			print(sp1[0], len(related_position))
			for rel_var1 in related_variant:
				rel_var1_string = chr1 + '-' + str(variant_position_dict[chr1][rel_var1])
				print(rel_var1_string)
				if mapping_to_var[rel_var1_string] != -1:
					matrix_id = mapping_to_var[rel_var1_string]
					if rel_var1_string in matrix[matrix_id][0]:
						row_id = matrix[matrix_id][0].index(rel_var1_string)
					else:
						continue
					matrix[matrix_id][1][row_id][row_id] += 1
				else:
					continue
				for rel_var2 in related_variant:
					rel_var2_string = chr1 + '-' + str(variant_position_dict[chr1][rel_var2])
					if mapping_to_var[rel_var2_string] != -1 and rel_var2_string != rel_var1_string:
						if rel_var2_string in matrix[matrix_id][0]:
							col_id = matrix[matrix_id][0].index(rel_var2_string)
						else:
							continue
						matrix[matrix_id][1][row_id][col_id] += 1
					else:
						continue
			for i_rel_var in range(len(related_variant)):
				i_var = related_variant[i_rel_var]
				VARIANT = [chr1, variant_position_dict[chr1][i_var], context_dict[chr1][i_var]]
				EXTRA = []
				homopolymer_on_ref = extract_homopolymer_on_reference(context_dict[chr1][i_var])
				EXTRA.append(homopolymer_on_ref)	
				quality1 = extract_base_quality(int(variant_position_dict[chr1][i_var]), int(pos1), qualstring1, CIGAR1)
				quality0 = extract_base_quality(int(variant_position_dict[chr1][i_var]), int(pos0), qualstring0, CIGAR0)
				base1 = extract_reference_base_read_base(int(variant_position_dict[chr1][i_var]), int(pos1), basestring1, context_dict[chr1][i_var], CIGAR1)
				base0 = extract_reference_base_read_base(int(variant_position_dict[chr1][i_var]), int(pos0), basestring0, context_dict[chr1][i_var], CIGAR0)
				cigar1 = extract_cigar_base(int(pos1), int(variant_position_dict[chr1][i_var]), CIGAR1, len(basestring1))
				cigar0 = extract_cigar_base(int(pos0), int(variant_position_dict[chr1][i_var]), CIGAR0, len(basestring0))
				indel1 = extract_nearby_indel(int(variant_position_dict[chr1][i_var]), int(pos1), CIGAR1)
				indel0 = extract_nearby_indel(int(variant_position_dict[chr1][i_var]), int(pos0), CIGAR0)
				homopolymer1 = extract_homopolymer_on_read(int(variant_position_dict[chr1][i_var]), int(pos1), basestring1, CIGAR1)
				homopolymer0 = extract_homopolymer_on_read(int(variant_position_dict[chr1][i_var]), int(pos0), basestring0, CIGAR0)
				INFOA_1 = INFO_1 + quality1 + base1 + cigar1 + [indel1, homopolymer1]
				INFOA_0 = INFO_0 + quality0 + base0 + cigar0 + [indel0, homopolymer0]
				if related_read[i_rel_var] == 0:
					RELATIONA = RELATION +  [strand0, strand1, flag0, flag1]
					out_list = INFOA_0 + RELATIONA + INFOA_1 + EXTRA + SUP_INFO_0 + SUP_INFO_1 + VARIANT + [related_position[i_rel_var], len(basestring0) - related_position[i_rel_var], CIGAR0.find('S')]
					out_list = [str(i) for i in out_list]
					if CIGAR0.find('S') > 0 or CIGAR0.find('H') > 0:
						continue
					OUTPUT.write('\t'.join(out_list) + '\n')	
				else:
					RELATIONA = RELATION + [strand1, strand0, flag1, flag0]
					out_list = INFOA_1 + RELATIONA + INFOA_0 + EXTRA + SUP_INFO_1 + SUP_INFO_0 + VARIANT + [related_position[i_rel_var], len(basestring1) - related_position[i_rel_var], CIGAR1.find('S')]
					out_list = [str(i) for i in out_list]
					if CIGAR1.find('S') > 0 or CIGAR1.find('H') > 0:
						continue
					OUTPUT.write('\t'.join(out_list) + '\n')
		else:
			size = -1
			samechr = 0
			RELATION = [str(size), str(samechr)]
			related_variant0 = []
			related_position0 = []
			related_variant1 = []
			related_position1 = []
			for i_var in range(len(variant_position_dict[chr1])):
				var_pos = variant_position_dict[chr1][i_var]
				var_loc_on_read1 = find_location_on_read(int(pos1), int(var_pos), CIGAR1, len(basestring1))
				if var_loc_on_read1 != "NA":
				#	print var_loc_on_read1
					if var_loc_on_read1 == "D":
#						related_variant1.append(i_var)
#						related_position1.append(var_loc_on_read1)
						continue
					elif basestring1[var_loc_on_read1] == variant_base_dict[chr1][i_var]:
						related_variant1.append(i_var)
						related_position1.append(var_loc_on_read1)
						read_pos_loc[var_pos].append(var_loc_on_read1)
						continue
			for rel_var1 in related_variant1:
				rel_var1_string = chr1 + '-' + str(variant_position_dict[chr1][rel_var1])
				if mapping_to_var[rel_var1_string] != -1:
					matrix_id = mapping_to_var[rel_var1_string]
					if rel_var1_string in matrix[matrix_id][0]:
						row_id = matrix[matrix_id][0].index(rel_var1_string)
					else:
						continue
					matrix[matrix_id][1][row_id][row_id] += 1
				else:
					continue
				for rel_var2 in related_variant1:
					rel_var2_string = chr1 + '-' + str(variant_position_dict[chr1][rel_var2])
					if mapping_to_var[rel_var2_string] != -1 and rel_var2_string != rel_var1_string:
						if rel_var2_string in matrix[matrix_id][0]:
							col_id = matrix[matrix_id][0].index(rel_var2_string)
						else:
							continue
						matrix[matrix_id][1][row_id][col_id] += 1
					else:
						continue
			for i_var in range(len(variant_position_dict[chr0])):
				var_pos = variant_position_dict[chr0][i_var]
				var_loc_on_read0 = find_location_on_read(int(pos0), int(var_pos), CIGAR0, len(basestring0))
				if var_loc_on_read0 != "NA":
				#	print var_loc_on_read0
					if var_loc_on_read0 == "D":
#						related_variant0.append(i_var)
#						related_position0.append(var_loc_on_read0)
						continue
					elif basestring0[var_loc_on_read0] == variant_base_dict[chr0][i_var]:
						related_variant0.append(i_var)
						related_position0.append(var_loc_on_read0)
						continue
			for rel_var1 in related_variant0:
				rel_var1_string = chr0 + '-' + str(variant_position_dict[chr0][rel_var1])
				if mapping_to_var[rel_var1_string] != -1:
					matrix_id = mapping_to_var[rel_var1_string]
					if rel_var1_string in matrix[matrix_id][0]:
						row_id = matrix[matrix_id][0].index(rel_var1_string)
					else:
						continue
					matrix[matrix_id][1][row_id][row_id] += 1
				else:
					continue
				for rel_var2 in related_variant0:
					rel_var2_string = chr0 + '-' + str(variant_position_dict[chr0][rel_var2])
					if mapping_to_var[rel_var2_string] != -1 and rel_var2_string != rel_var1_string:
						if rel_var2_string in matrix[matrix_id][0]:
							col_id = matrix[matrix_id][0].index(rel_var2_string)
						else:
							continue
						matrix[matrix_id][1][row_id][col_id] += 1
					else:
						continue
			print(sp1[0], len(related_variant1), len(related_variant0))
			for i_rel_var in range(len(related_variant1)):
				i_var = related_variant1[i_rel_var]
				VARIANT = [chr1, variant_position_dict[chr1][i_var], context_dict[chr1][i_var]]
				EXTRA = []
				homopolymer_on_ref = extract_homopolymer_on_reference(context_dict[chr1][i_var])
				EXTRA.append(homopolymer_on_ref)
				quality1 = extract_base_quality(int(variant_position_dict[chr1][i_var]), int(pos1), qualstring1, CIGAR1)
				quality0 = extract_base_quality(int(variant_position_dict[chr1][i_var]), int(pos0), qualstring0, CIGAR0)
				base1 = extract_reference_base_read_base(int(variant_position_dict[chr1][i_var]), int(pos1), basestring1, context_dict[chr1][i_var], CIGAR1)
				base0 = extract_reference_base_read_base(int(variant_position_dict[chr1][i_var]), int(pos0), basestring0, context_dict[chr1][i_var], CIGAR0)
				cigar1 = extract_cigar_base(int(pos1), int(variant_position_dict[chr1][i_var]), CIGAR1, len(basestring1))
				cigar0 = extract_cigar_base(int(pos0), int(variant_position_dict[chr1][i_var]), CIGAR0, len(basestring0))
				indel1 = extract_nearby_indel(int(variant_position_dict[chr1][i_var]), int(pos1), CIGAR1)
				indel0 = extract_nearby_indel(int(variant_position_dict[chr1][i_var]), int(pos0), CIGAR0)
				homopolymer1 = extract_homopolymer_on_read(int(variant_position_dict[chr1][i_var]), int(pos1), basestring1, CIGAR1)
				homopolymer0 = extract_homopolymer_on_read(int(variant_position_dict[chr1][i_var]), int(pos0), basestring0, CIGAR0)
				INFOA_1 = INFO_1 + quality1 + base1 + cigar1 + [indel1, homopolymer1]
				INFOA_0 = INFO_0 + quality0 + base0 + cigar0 + [indel0, homopolymer0]
				RELATIONA = RELATION + [strand1, strand0, flag1, flag0]
				out_list = INFOA_1 + RELATIONA + INFOA_0 + EXTRA + SUP_INFO_1 + SUP_INFO_0 + VARIANT + [related_position1[i_rel_var], len(basestring1) - related_position1[i_rel_var], CIGAR1.find('S')]
				out_list = [str(i) for i in out_list]
				if CIGAR1.find('S') > 0 or CIGAR1.find('H') > 0:
					continue
				OUTPUT.write('\t'.join(out_list) + '\n')
			for i_rel_var in range(len(related_variant0)):
				i_var = related_variant0[i_rel_var]
				VARIANT = [chr0, variant_position_dict[chr0][i_var], context_dict[chr0][i_var]]
				EXTRA = []
				homopolymer_on_ref = extract_homopolymer_on_reference(context_dict[chr0][i_var])
				EXTRA.append(homopolymer_on_ref)
				quality1 = extract_base_quality(int(variant_position_dict[chr0][i_var]), int(pos1), qualstring1, CIGAR1)
				quality0 = extract_base_quality(int(variant_position_dict[chr0][i_var]), int(pos0), qualstring0, CIGAR0)
				base1 = extract_reference_base_read_base(int(variant_position_dict[chr0][i_var]), int(pos1), basestring1, context_dict[chr0][i_var], CIGAR1)
				base0 = extract_reference_base_read_base(int(variant_position_dict[chr0][i_var]), int(pos0), basestring0, context_dict[chr0][i_var], CIGAR0)
				cigar1 = extract_cigar_base(int(pos1), int(variant_position_dict[chr0][i_var]), CIGAR1, len(basestring1))
				cigar0 = extract_cigar_base(int(pos0), int(variant_position_dict[chr0][i_var]), CIGAR0, len(basestring0))
				indel1 = extract_nearby_indel(int(variant_position_dict[chr0][i_var]), int(pos1), CIGAR1)
				indel0 = extract_nearby_indel(int(variant_position_dict[chr0][i_var]), int(pos0), CIGAR0)
				homopolymer1 = extract_homopolymer_on_read(int(variant_position_dict[chr0][i_var]), int(pos1), basestring1, CIGAR1)
				homopolymer0 = extract_homopolymer_on_read(int(variant_position_dict[chr0][i_var]), int(pos0), basestring0, CIGAR0)
				INFOA_1 = INFO_1 + quality1 + base1 + cigar1 + [indel1, homopolymer1]
				INFOA_0 = INFO_0 + quality0 + base0 + cigar0 + [indel0, homopolymer0]
				RELATIONA = RELATION + [strand0, strand1, flag0, flag1]
				out_list = INFOA_0 + RELATIONA + INFOA_1 + EXTRA + SUP_INFO_0 + SUP_INFO_1 + VARIANT + [related_position0[i_rel_var], len(basestring0) - related_position0[i_rel_var], CIGAR0.find('S')]
				out_list = [str(i) for i in out_list]
				if CIGAR0.find('S') > 0 or CIGAR0.find('H') > 0:
					continue
				OUTPUT.write('\t'.join(out_list) + '\n')


OUTPUT.close()


# In[36]:

cluster_var = open(sys.argv[4] + ".var.txt", 'w')
CLUSTERED_FRAC = 0.9
NONCLUSTERED_COUNT = 2
clustered_variants = []
for item in matrix:
	var_list = item[0]
	cluster_var.write('\t'.join(var_list) + '\n')
	co_matrix = item[1]
	for i_var_main in range(len(var_list)):
		cnt_var_main = co_matrix[i_var_main][i_var_main]
		if cnt_var_main == 0:
			print(var_list[i_var_main])
			continue
		flag = 0
		pos = int(var_list[i_var_main].split('-')[1])
		for i_var_else in range(len(var_list)):
			if i_var_main == i_var_else:
				continue
			pos_else = int(var_list[i_var_else].split('-')[1])
			if abs(pos_else - pos) < 3:
				continue
			cnt_cooccur = co_matrix[i_var_main][i_var_else]
			if float(cnt_cooccur)/float(cnt_var_main) > CLUSTERED_FRAC and (cnt_var_main - cnt_cooccur) < NONCLUSTERED_COUNT:
				flag = 1
		if flag == 1:
			clustered_variants.append(var_list[i_var_main])

cluster_var.close()


# In[39]:


clusters = open(sys.argv[4], 'w')
for var in clustered_variants:
	sp = var.split('-')
	sp.append(sp[1])
	sp[1] = str(int(sp[1])-1)
	clusters.write('\t'.join(sp) + '\n')


variant_position_chrom = variant_position_dict.keys()
for chrom in variant_position_chrom:
	for var_pos in variant_position_dict[chrom]:
		if len(set(read_pos_loc[var_pos])) < 0:
			sp = [chrom, str(var_pos)]
			sp.append(sp[1])
			sp[1] = str(int(sp[1])-1)
			clusters.write('\t'.join(sp) + '\n')
		if len(set(read_pos_loc[var_pos])) == 0:
			print(chrom, var_pos)
clusters.close()

