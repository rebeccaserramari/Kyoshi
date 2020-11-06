configfile: "config.yaml"


chromosomes = config["chromosomes"]
samples = config["samples"]
haplos = config["haplos"]
pops = config["haplos"]
callset_file = config["callset_file"]
callset_file_unzipped = config["callset_file_unzipped"]
callset_folder = config["callset_folder"]
results = config["results"]
panel = config["panel"]
populations = config["populationfile"]
result_plots = config["result_plots"]
num_of_samples = config["NSAMPLES"]
bootstrap = config["BOOTSTRAP"]


rule all:
	input:
		expand(result_plots + "{sample}/{chromosome}/{sample}_{haplotype}_{chromosome}.pdf",chromosome=chromosomes, sample=samples,haplotype=haplos)	

rule plot:
	input:
		rules.average_probabilities.output
	output:
		result_plots + "{sample}/{chromosome}/{sample}_{haplotype}_{chromosome}.pdf"
	shell:
		"Rscript plotting.r {input} {output}"

rule subset_vcf:
	input:
		callset_folder + callset_file
	output:
		callset_folder + "{sample}/{sample}." + callset_file
	threads:
		40
	shell:
		"bcftools view --threads {threads} {input} -s {wildcards.sample} -O z -o {output}"

rule index_vcf:
	input:
		rules.subset_vcf.output
	output:
		callset_folder + "{sample}/{sample}" + callset_file + ".csi"
	shell:
		"bcftools index {input}"

rule subset_sample:
	input:
		vcf = rules.subset_vcf.output,
		index = rules.index_vcf.output
	output:
		callset_folder + "{sample}/{sample}.{chromosome}." + callset_file_unzipped
	threads:
		40
	shell:
		"bcftools view --threads {threads} {input.vcf} -r chr{wildcards.chromosome} -O v -o {output}"

rule sample:
	input:
		pop = populations
	output:
		samples = results+"{sample}/samples/samples_bs_{bs}"
	params:
		number = num_of_samples
	shell:
		"/usr/bin/time -v python scripts/sampling.py"		
		
rule create_populationfile:
	input:
		pop = populations,
		samples = results + "{sample}/samples/samples_bs_{bs}"
	output:
		population = results + "{sample}/populations/population_{bs}.tsv"
	shell:
		"/usr/bin/time -v python scripts/write_populationfile.py"
			
rule compute_fwdbkw:
	input:
		vcf = rules.subset_sample.output,
		panel = panel,
		population = results+"{sample}/populations/population_{bs}.tsv"
	output:
		variants = results+"{sample}/{chromosome}/variants_{sample}_{chromosome}_{bs}.txt",
		posterior_h1 = results+"{sample}/{chromosome}/posterior_h1_{chromosome}_{bs}.txt",
		posterior_h2 = results+"{sample}/{chromosome}/posterior_h2_{chromosome}_{bs}.txt"
	shell:
		'/usr/bin/time -v python hapcolor.py {input.vcf} {input.panel} {input.population} {output.variants} {output.posterior_h1} {output.posterior_h2}'


rule postprocess:
	input:
		population = rules.create_populationfile.output,
		variants = results+"{sample}/{chromosome}/variants_{sample}_{chromosome}_{bs}.txt",
		posterior = results+"{sample}/{chromosome}/posterior_{haplo}_{chromosome}_{bs}.txt",
	output:
		prob = results+"{sample}/{chromosome}/{haplo}_probabilities_chr{chromosome}_{bs}.tsv"
		
	shell:
		'/usr/bin/time -v python scripts/postprocess.py'


rule average_probabilities:
	input:
		expand(results+"{{sample}}/{{chromosome}}/{{haplo}}_probabilities_chr{{chromosome}}_{bs}.tsv", bs=BOOTSTRAP)
	output:
		results+"{sample}/{chromosome}/{haplo}_averaged_probabilities_chr{chromosome}.tsv"
	run:
		variants = []
		with open(input[0]) as firstfile:
			for i,line in enumerate(firstfile):
				variants.append(line.strip().split('\t')[:2])
								
		amr, eur, afr, eas, sas = [0 for i in range(len(variants))],[0 for i in range(len(variants))],[0 for i in range(len(variants))],[0 for i in range(len(variants))],[0 for i in range(len(variants))]
		for filename in input:
			with open(filename) as infile:
				for i,line in enumerate(infile):
					afr[i] += float(line.strip().split('\t')[2])
					amr[i] += float(line.strip().split('\t')[3])
					eas[i] += float(line.strip().split('\t')[4])
					eur[i] += float(line.strip().split('\t')[5])
					sas[i] += float(line.strip().split('\t')[6])
		afr = [round(i/len(BOOTSTRAP),4) for i in afr]
		amr = [round(i/len(BOOTSTRAP),4) for i in amr]
		eas = [round(i/len(BOOTSTRAP),4) for i in eas]
		eur = [round(i/len(BOOTSTRAP),4) for i in eur]
		sas = [round(i/len(BOOTSTRAP),4) for i in sas]
		with open(output[0],'w') as out:
			for i in range(len(variants)):
				out.write('\t'.join(variants[i])+'\t'+str(afr[i])+'\t'+str(amr[i])+'\t'+str(eas[i])+'\t'+str(eur[i])+'\t'+str(sas[i])+'\n')
				