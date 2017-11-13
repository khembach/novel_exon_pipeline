
##################

# RSEM

##################
# RSEM version 1.3.0

rule prepare_reference:
	input:
		"{GENOMEDIR}*.fa",
		GTF
	output: "{RSEMREF}.n2g.idx.fa"
	threads: CORES
	shell:
		"rsem-prepare-reference --gtf {GTF} --star --star-sjdboverhang 100 -p {CORES} {GENOMEDIR} {REFDIR}"
## we have paired-end stranded reads from Illumina HiSeq 2000 (stranded TruSeq libary preparation with dUTPs )
## --> first read comes from the reverse strand: set --strandedness reverse

rule calculate_expression:
	input:
		fa1 = FASTQ1,
		fa2 = FASTQ2
	output: "simulation/{SAMPLENAME}.isoforms.results"
	threads: CORES
	shell:
		"rsem-calculate-expression --paired-end {input.fa1} {input.fa2} {RSEMREF} simulation/{SAMPLENAME} --strandedness reverse -p {CORES} --star --star-gzipped-read-file"
		#--fragment-lenth-min 200   --> should we set this parameter?


### modify the RSEM model file (remove the probability of generating reads with quality 2 to prevent sampling of reads with overall low quality)
rule parse_markov_prob:
	input: "simulation/{SAMPLENAME}.stat/{SAMPLENAME}.model"
	output: "simulation/{SAMPLENAME}.stat/{SAMPLENAME}_markov_prob"
	shell: "sed -n '14,114p;115q' {input} > {output}"

rule modify_markov_prob:
	input: "simulation/{SAMPLENAME}.stat/{SAMPLENAME}_markov_prob"
	output: "simulation/{SAMPLENAME}.stat/{SAMPLENAME}_markov_prob_modified"
	script: "scripts/modify_RSEM_quality_score_distribution.R"

rule modify_RSEM_model:
	input:
		model="simulation/{SAMPLENAME}.stat/{SAMPLENAME}.model",
		modified_markov_prob = "simulation/{SAMPLENAME}.stat/{SAMPLENAME}_markov_prob_modified"
	output: "simulation/{SAMPLENAME}.stat/{SAMPLENAME}.model_modified"
	shell:
		"sed -n '1,13p;14q' {input.model} > {output} \n"
		"cat {input.modified_markov_prob} >> {output} \n"
		"sed -n '115,$p' {input} >> {output}"

rule simulate_data:
	input: model = "simulation/{SAMPLENAME}.stat/{SAMPLENAME}.model_modified".format(SAMPLENAME=SAMPLENAME),
		iso = "simulation/{SAMPLENAME}.isoforms.results".format(SAMPLENAME=SAMPLENAME)
	output: protected( expand("simulation/simulated_data/simulated_reads_chr19_22_{nr}.fq", nr = [1,2]) )
	shell:
		"rsem-simulate-reads {RSEMREF} {input.model} {input.iso} {THETA} {N_READS} simulation/simulated_data/simulated_reads_chr19_22 --seed {SEED}"
