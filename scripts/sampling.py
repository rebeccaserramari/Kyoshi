import sys
import math
import random
from collections import defaultdict
from itertools import chain, groupby
import statistics
import pdb

sampling([snakemake.input["pop"], snakemake.params["number"], snakemake.output["samples"])
	
def sampling(popfile, num_samples, outputfile):
	mapping = defaultdict(list)
	supermapping = defaultdict(list)
	allsamples = []
	with open(popfile) as f:
		for i,line in enumerate(f):
			if (i>0):
				sample = line.strip().split('\t')[0]
				superpop = line.strip().split('\t')[2]
				
				supermapping[superpop].append(sample)

	total_samples = []
	#sample num_samples from each superpopulation present in the population file
	for super_pop in supermapping.keys():
		tmp_sampling = random.choices(supermapping[super_pop], k = num_samples)
		total_samples.extend(tmp_sampling)	
	
	assert(len(total_samples) == len(supermapping.keys())*num_samples)

	with open(outputfile, 'w') as outfile:
		outfile.write("\n".join(total_samples))
	return(total_samples)