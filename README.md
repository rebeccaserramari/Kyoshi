# Kyoshi
reference-based ancestry inference and colouring of haplotypes

author: Rebecca Serra Mari

## Description
HMM-based method to infer the local ancestry of each genomic region of the haplotype-resolved assemblies.
Takes an assembled haplotype and a phased reference panel as input and outputs the most likely ancestral population for each variant position.
The computation of the most likely ancestral population for each genomic position is performed by setting up a Hidden Markov Model that models the reference haplotypes from the panel as hidden state sequences and the target haplotype as the observed sequence and running a forward-backward algorithm on this HMM. 
At each variant position, costs are inferred for a) discrepancies between target and reference haplotype and b) switches to a different reference haplotype happening between two variant positions. The resulting probability will be high for reference haplotypes that largely coincide with the target, for blocks as long as possible. 

This computation results in a set of probabilities for each variant position, consisting of one probability per reference haplotype. The probabilities of reference samples that belong to the same superpopulation are added to get a likelihood for every superpopulation at each position.

## Usage
Clone the repository into a separate conda environment and run `snakemake`. Adapt the file `config.yaml` with the correct input file paths.
   
