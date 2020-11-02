import sys
import logging
import itertools
import operator
from cyvcf2 import VCF, Writer
from collections import defaultdict
import collections
import gc
import math
from sys import getsizeof, stderr
import numpy as np
from difflib import SequenceMatcher


	
#reads a VCF file and extracts the haplotype information, both haplotypes are output as lists of ints
def compute_haplotypes(var_map, variant_set):
	haplotypes = ["",""]

	for index in var_map.keys():
		v = var_map[index]
		counter = 0
		for genotype in v.genotypes:
		#	if (genotype[0]==0 and genotype[1]==0):
		#		haplotypes[counter] += '-'
		#		haplotypes[counter+1] += '-'
			if (genotype[0] == -1 or genotype[1] == -1):
				haplotypes[counter] += '-'
				haplotypes[counter+1] += '-'
			else:
				haplotypes[counter] += str(genotype[0])
				haplotypes[counter+1] += str(genotype[1])
			counter += 2
	assert(len(haplotypes) == 2)
	assert(len(haplotypes[0]) == len(var_map))

	hap1 = [int(i) for i in haplotypes[0]]
	hap2 = [int(i) for i in haplotypes[1]]

	return(hap1, hap2)

	
def compute_referencepanel(ref_file,  samples_to_use, variants_to_use):
	#compute size of ref_matrix and initialize it with zeros
	#width=number of variants between first and last relevant variant
	#height=twice the number of samples, as each sample offers two haplotypes

	height = 2*len(samples_to_use)

	refmatrix = ['' for i in range(height)]
	vcf_ref = VCF(ref_file, lazy=True)
	doubleset = set()
	index = -1
	homo = 0
	homozyg = False
	for variant in variants_to_use:
		if (variant.POS not in doubleset):
			index += 1
			doubleset.add(variant.POS)
			counter = 0
			homozyg = False
			for baseindex in samples_to_use:		
				genotype = variant.genotypes[baseindex]
				if (genotype[0]==-1 or genotype[1]==-1):
					refmatrix[counter] += '-'
					refmatrix[counter+1] += '-'
			#	elif (genotype[0]==0 and genotype[1]==0):
			#		refmatrix[counter] += '-'
			#		refmatrix[counter+1] += '-'
				else:
					if (variant.gt_phases[baseindex] == False and variant.gt_types[baseindex] == 1):
						refmatrix[counter] += '-'
						refmatrix[counter+1] += '-'
					else:					
						refmatrix[counter] += str(genotype[0])
						refmatrix[counter+1] += str(genotype[1])
				if (genotype[0] == genotype[1]):
					homozyg = True
				counter += 2
			if homozyg:
				homo += 1
	return (refmatrix)
	

