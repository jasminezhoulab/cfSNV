from parameter import *
import numpy as np

def string_to_base(basestring):
# input: a base string
# output: list of bases in the base string
# string type in python is iterable, so output is the same as input here
	return basestring

def string_to_qual(qualitystring):
# input: a qual string
# output: a list of error probabilities
# qual is the phred score of error probability
	prob = np.array([10.0**(-float(ord(i)-33)/10.0) for i in qualitystring])
	return prob


def filter_low_qual(basestring, quallist, maplist):
	id = ( quallist < BASEQUAL_THRESHOLD ) * ( maplist < MAPQUAL_THRESHOLD )
	new_basestring = "".join([basestring[i] for i in range(len(basestring)) if id[i] == True])
	new_quallist = quallist[id]
	new_maplist = maplist[id]
	return new_basestring, new_quallist, new_maplist

def filter_low_qual_with_string(basestring, quallist, maplist, qualstring, mapstring):
	id = ( quallist < BASEQUAL_THRESHOLD ) * ( maplist < MAPQUAL_THRESHOLD )
	new_basestring = "".join([basestring[i] for i in range(len(basestring)) if id[i] == True])
	new_quallist = quallist[id]
	new_maplist = maplist[id]
	new_qualstring = "".join([qualstring[i] for i in range(len(qualstring)) if id[i] == True])
	new_mapstring = "".join([mapstring[i] for i in range(len(mapstring)) if id[i] == True])
	return new_basestring, new_quallist, new_maplist, new_qualstring, new_mapstring



def find_major_variant(basecount_nc, basecount_ef):
# input: converted basestring
# output: variant with the max count of observations
#       haven't thought of the case where the max count corresponds to a variant with strong strand bias
	alternative_nucleotide = ['A', 'C', 'T', 'G']
	count = []
	for nucleotide in alternative_nucleotide:
		count.append( basecount_nc[nucleotide] + basecount_nc[nucleotide.lower()] + 2*basecount_ef[nucleotide] + 2*basecount_ef[nucleotide.lower()] )
	max_count = max(count)
	variant_base = alternative_nucleotide[count.index(max_count)]
	return variant_base


def observe_variant_probability(tumor_fraction, joint_genotype):
# input: tumor fraction, joint genotype
# output: theoretical allele frequencies of reference allele and non-reference allele
# combine tumor fraction and joint genotype to the frequencies of nucleotides
	variant_allele_observed_probability = dict([])
	if joint_genotype == 'AA/AA':
		variant_allele_observed_probability['R'] = 1.0
		variant_allele_observed_probability['X'] = 0.0
	elif joint_genotype == 'AA/AB':
		variant_allele_observed_probability['R'] = 1.0 - 0.5*tumor_fraction
		variant_allele_observed_probability['X'] = 0.5*tumor_fraction
	elif joint_genotype == 'AA/BB':
		variant_allele_observed_probability['R'] = 1.0 - tumor_fraction
		variant_allele_observed_probability['X'] = tumor_fraction
	elif joint_genotype == 'AB/AA':
		variant_allele_observed_probability['R'] = 0.5 + 0.5*tumor_fraction
		variant_allele_observed_probability['X'] = 0.5 - 0.5*tumor_fraction
	elif joint_genotype == 'AB/AB':
		variant_allele_observed_probability['R'] = 0.5
		variant_allele_observed_probability['X'] = 0.5
	elif joint_genotype == 'AB/BB':
		variant_allele_observed_probability['R'] = 0.5 - 0.5*tumor_fraction
		variant_allele_observed_probability['X'] = 0.5 + 0.5*tumor_fraction
	elif joint_genotype == 'BB/AA':
		variant_allele_observed_probability['R'] = tumor_fraction
		variant_allele_observed_probability['X'] = 1.0 - tumor_fraction
	elif joint_genotype == 'BB/AB':
		variant_allele_observed_probability['R'] = 0.5*tumor_fraction
		variant_allele_observed_probability['X'] = 1.0 - 0.5*tumor_fraction
	elif joint_genotype == 'BB/BB':
		variant_allele_observed_probability['R'] = 0.0
		variant_allele_observed_probability['X'] = 1.0
	else:
		print("error with input genotype")
	return variant_allele_observed_probability


def translate_basestring(basestring, variant_base):
# input: basestring in rebase file
# output: non-reference alleles, converted basestring, base count string
# find all observed non-reference alleles
# count occurrence of non-reference alleles
# convert reference alleles to 'R', and non-reference alleles to 'X'
	newstring = basestring
	newstring = newstring.replace(variant_base, 'X')
	newstring = newstring.replace(variant_base.lower(), 'X')
	newstring = newstring.replace('r', 'R')
	return np.array(list(newstring), dtype = str)


def count_base(basestring):
# input: basestring in rebase file
# output: base count dict
# count occurrence of non-reference alleles
	basecount = {}
	for i in "AaCcGgTtRr":
		basecount[i] = basestring.count(i)
	return basecount



def convert_to_Decimal(x):
# input: integer or float
# output: decimal object
	return Decimal(str(x))



def calculate_joint_genotype_tumor_fraction_loglikelihood(tumor_fraction, basestring, quallist, maplist, joint_genotype, variant_base):
# input: tumor fraction, converted basestring, base quality list, joint genotype
# output: log likelihood log( P( X | theta, joint_genotype ) )
	newbasestring = translate_basestring(basestring, variant_base)
	interest_id = (newbasestring == 'R') + (newbasestring == 'X')
	VAF = observe_variant_probability(tumor_fraction, joint_genotype)
	#print VAF, newbasestring, basestring, variant_base,  interest_id
#	print len(basestring), len(newbasestring), len(quallist), len(maplist)
	interest_quallist = quallist[interest_id]
	interest_maplist = maplist[interest_id]
	interest_VAFlist = np.array([VAF[i] for i in newbasestring[interest_id]])
	#print interest_VAFlist, interest_quallist, interest_maplist 
	loglikelihood = sum(np.log(interest_VAFlist*(1-interest_quallist)*(1-interest_maplist) + (1-interest_VAFlist)*(interest_quallist+interest_maplist+interest_quallist*interest_maplist)))
	return loglikelihood



