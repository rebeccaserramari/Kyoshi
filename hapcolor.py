"""
Colouring chromosomes by ancestries inferred by a reference panel

"""
from __future__ import absolute_import
import logging
from utils import compute_referencepanel, compute_haplotypes
import sys
import argparse
from cyvcf2 import VCF, Writer
import platform
import resource
import random
from core import hmm
from collections import defaultdict
from operator import itemgetter


logger = logging.getLogger(__name__)


	
def validate(args, parser):
	pass

def genotype_exists(v):
	for genotype in v.genotypes:
		if (genotype[0] == -1 or genotype[1] == -1):
			return(False)
	return(True)
	
def included_variants(ref_file, variant_file):
	ref_map = defaultdict()
	chromosome_reference_set = set()
	refvcf = VCF(ref_file)
	for variant in refvcf:
		ref_map[variant.POS] = variant
		chromosome_reference_set.add(variant.CHROM)
	ref_set = set(ref_map.keys())

	#for multiple chromosomes, raise an error
	if (len(chromosome_reference_set)>1):
		logger.error("The given reference file is supposed to contain variants from one chromosome only!")
		sys.exit(1)
	#create the list of variants to be considered	
	chromosome_variant_set = set()
	vcf = VCF(variant_file)
	doubleset = set()
	var_map = {}
	for v in vcf:
		#use variant if it is present in the panel, heterozygous and phased	
		if (v.gt_types[0] == 1) and (v.gt_phases[0] == True) and (v.is_snp) and (v.POS not in doubleset) and genotype_exists(v):	
			doubleset.add(v.POS)
			if (v.CHROM not in chromosome_reference_set and ("chr"+v.CHROM) not in chromosome_reference_set and v.CHROM.replace('chr','') not in chromosome_reference_set):
				logger.error("Reference and target file must contain variants from the same chromosome!")
				sys.exit(1)
			chromosome_variant_set.add(v.CHROM)
			if (v.POS in ref_set):
				var_map[v.POS] = v
			
	variant_set = set(var_map.keys())

	if (len(chromosome_variant_set)>1):
		logger.error("The given variant file is supposed to contain variants from one chromosome only!")
		sys.exit(1)
	if not variant_set:
		logger.error("No positions found in both the reference and the target file.")
		sys.exit(1)
		
	return(variant_set, ref_map, var_map)


			
def run_connect(variant_file, reference_file, population_file, output_variants, file_posterior_h1, file_posterior_h2, mapfile=None):
	
	logger.info("Choosing suitable variant positions.")
	(variant_set, vcf_ref_map, var_map) = included_variants(reference_file, variant_file)


	variants_to_use = [value for key,value in vcf_ref_map.items() if key in variant_set]
	positionslist = [variants_to_use[var].POS for var in range(len(variants_to_use))]
	positionslist.sort()

	with open(output_variants,'w') as variants:
		for pos in range(len(positionslist)):
			if (pos == len(positionslist)-1):
				variants.write(str(positionslist[pos]))
			else:
				variants.write(str(positionslist[pos])+',')
	variants.close()

	recomb_map = {}
	recomb_list = []
	if mapfile:
		genetic_map = defaultdict()	
		with open(mapfile) as f:
			for line in f:
				if line.startswith('pos'):
					continue
				else: 
					position = line.strip().split('\t')[0]
					cm = line.strip().split('\t')[2]
					genetic_map[int(position)] = float(cm)

		in_map, not_inmap = 0, 0
		final_map = {}
		vcount = 0
		for i in range(0, len(variants_to_use)):
			v = variants_to_use[i].POS
			if (v in set(genetic_map.keys())):
				in_map += 1
				final_map[v] = genetic_map[v]
			else:
				prev_elements = sorted([j for j in sorted(genetic_map.keys()) if j < v], reverse=True)
				succ_elements = sorted([j for j in sorted(genetic_map.keys()) if j > v])
				not_inmap += 1	
				if (len(prev_elements) == 0):
					succ = succ_elements[1]
					pred = succ_elements[0]
				elif (len(succ_elements) == 0):
					succ = prev_elements[0]
					pred = prev_elements[1]
				else:
					succ = succ_elements[0]
					pred = prev_elements[0] 
				diff_pos = succ-pred
				assert(diff_pos > 0)
				diff_cm = genetic_map[succ]-genetic_map[pred]
				assert(diff_cm >= 0)
				current_diff = abs(v - pred)
				assert(current_diff >= 0)
				interpolated_cm = (diff_cm/diff_pos)*current_diff + genetic_map[pred]
				assert(interpolated_cm >= 0)
				final_map[v] = interpolated_cm
		vcount = 0
		for i in range(len(variants_to_use)):
			v = variants_to_use[i].POS
			if (i == 0):
				rf = final_map[v]
			else:
				u = variants_to_use[i-1].POS
				rf = final_map[v] - final_map[u]
			
			recomb_map[v] = rf
			recomb_list.append(rf)
	else:
		for i in range(len(variants_to_use)):
			v = variants_to_use[i].POS
			recomb_map[v] = 0
			recomb_list.append(0)
			
	#writes the original haplotypes from the target into strings
	(haplo1,haplo2) = compute_haplotypes(var_map, variant_set)
	
	vcf_ref = VCF(reference_file, lazy=True)
	samples_to_use = list(range(len(vcf_ref.samples)))
	refmatrix = []

	samples = []
	with open(population_file) as pops:
		for i,line in enumerate(pops):
			samples.append(line.strip().split('\t')[0])
	
	indices = []		
	for sam in samples:
		if sam in vcf_ref.samples:
			indices.append(vcf_ref.samples.index(sam))

	logger.info("Computing reference panel")
	refmatrix = compute_referencepanel(reference_file, indices, variants_to_use)
		
	E_whole = []

	#performs the forward backward computation on the HMM
	logger.info("Dynamic programming: Performing haplotype block parsing")
	(path1, path2) = hmm(haplo1, haplo2, E_whole, refmatrix, 1,  recomb_list, file_posterior_h1, file_posterior_h2)

	logger.info('\n== SUMMARY ==')
	if sys.platform == 'linux':
		memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		logger.info('Maximum memory usage: %.3f GB', memory_kb / 1E6)

def main(args):
	run_connect(**vars(args))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('variant_file', metavar='BLOCKS', help='VCF file of target sample')
	parser.add_argument('reference_file', metavar='REFERENCE', help='VCF file that serves as reference panel, this will be used to infer ancestries')
	parser.add_argument('population_file', metavar='SAMPLES', help='File containing the samples and population information. ')
	parser.add_argument('output_variants', metavar='OUTPUT', help='Writes the considered variants into output file for later use for plotting.')
	parser.add_argument('file_posterior_h1', metavar='POSTERIOR_H1', help='Output file storing the posterior probabilities for haplotype 1.')
	parser.add_argument('file_posterior_h2', metavar='POSTERIOR_H2', help='Output file storing the posterior probabilities for haplotype 2.')
	parser.add_argument('--mapfile', metavar="MAP", help='Genetic map')
	args = parser.parse_args()
	
	main(args)		