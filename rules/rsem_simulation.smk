
##################

# RSEM

##################
# RSEM version 1.3.0

rule prepare_reference:
	input:
		gtf = config["gtf"]
	output: 
		config["RSEMREF"] + ".n2g.idx.fa"
	params:
		gdir = config["GENOMEDIR"],
	threads: 
		config["cores"]
	shell:
		"rsem-prepare-reference --gtf {input.gtf} --star --star-sjdboverhang 100 -p {threads} {params.gdir} {output}"
## we have paired-end stranded reads from Illumina HiSeq 2000 (stranded TruSeq libary preparation with dUTPs )
## --> first read comes from the reverse strand: set --strandedness reverse

rule calculate_expression:
	input:
		fa1 = config["FASTQ1"],
		fa2 = config["FASTQ2"]
	output: 
		"simulation/" + config["SAMPLENAME"] + ".isoforms.results"
		"simulation/" + config["SAMPLENAME"] + ".stat/" + config["SAMPLENAME"] + ".model"
	threads: 
		config["cores"]
	params:
		samplename = config["SAMPLENAME"],
		rsemref = config["RSEMREF"]
	shell:
		"rsem-calculate-expression --paired-end {input.fa1} {input.fa2} {params.rsemref} simulation/{params.samplename} --strandedness reverse -p {threads} --star --star-gzipped-read-file"
		#--fragment-lenth-min 200   --> should we set this parameter?


### modify the RSEM model file (remove the probability of generating reads with quality 2 to prevent sampling of reads with overall low quality)
rule parse_markov_prob:
	input: 
		"simulation/{samplename}.stat/{samplename}.model"
	output: 
		"simulation/{samplename}.stat/{samplename}_markov_prob"
	shell: 
		"sed -n '14,114p;115q' {input} > {output}"

rule modify_markov_prob:
	input: 
		"simulation/{samplename}.stat/{samplename}_markov_prob"
	output: 
		"simulation/{samplename}.stat/{samplename}_markov_prob_modified"
	script: 
		"scripts/modify_RSEM_quality_score_distribution.R"

rule modify_RSEM_model:
	input:
		model = "simulation/{samplename}.stat/{samplename}.model",
		modified_markov_prob = "simulation/{samplename}.stat/{samplename}_markov_prob_modified"
	output: 
		"simulation/{samplename}.stat/{samplename}.model_modified"
	shell:
		"sed -n '1,13p;14q' {input.model} > {output} \n"
		"cat {input.modified_markov_prob} >> {output} \n"
		"sed -n '115,$p' {input} >> {output}"

rule simulate_data:
	input: 
		model = "simulation/" + config["SAMPLENAME"] + ".stat/" + config["SAMPLENAME"] + ".model_modified",
		iso = "simulation/" + config["SAMPLENAME"] + ".isoforms.results"
	output: 
		protected("simulation/simulated_data/simulated_reads_chr19_22_1.fq"),
		protected("simulation/simulated_data/simulated_reads_chr19_22_2.fq"),
	params:
		rsemref = config["RSEMREF"],
		theta = config["THETA"],
		n_reads = config["N_READS"],
		seed = config["SEED"]
	shell:
		"rsem-simulate-reads {params.rsemref} {input.model} {input.iso} {params.theta} {params.n_reads} simulation/simulated_data/simulated_reads_chr19_22 --seed {params.seed}"
