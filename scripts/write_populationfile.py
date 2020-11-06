import sys
import math
import random
from collections import defaultdict


def write_populationfile(samplefile, popfile, outfile):
	with open(samplefile) as sams:
		samplenames = sams.read().splitlines()
	lines = []
	with open(popfile) as f:
		for i,line in enumerate(f):
			if (line.strip().split()[0] in samplenames):
				lines.append(line)
	
	with open(outfile,'w') as ofile:
		for line in lines:
			ofile.write(line)
	return()


if __name__ == '__main__':
	write_populationfile(snakemake.input["samples"], snakemake.input["pop"], snakemake.output["population"])	
	