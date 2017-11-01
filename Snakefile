#######
## The directory contains a .Renviron file that specifies the R library. This is necessary, because Rscript does not find the libraries without.
## Is there a way to specify the R, Rscript version and the library in the snakemake file?



## data was downloaded from ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX160/SRX1603411/SRR3192428/
## GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2072377
FASTQ1 = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/SRR3192428_1.fastq.gz"
FASTQ2 = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/SRR3192428_2.fastq.gz"
GTF = "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
CORES = 25
GENOMEDIR = "/home/Shared/data/seq/microexon_simulation/microexon_study/genome/reference/"
RSEMREF = "reference/RSEM/hg19_chr19_22"
REALSAMPLE = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/"
SIMDIR = "simulation/"
SAMPLENAME = "SRR3192428"
THETA = 0.05
N_READS = 1118017
## the number of reads that were mapped to chromosome 19 and 22: 1118017
## 8306564 = number of reads = [348875684]/4*{1500+500}/(21000)        [lines in real fastq data] {genes in chr19 + chr22} / (total genes)
SEED = 584


rule all:
	input:
		"{SIMDIR}reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf".format(SIMDIR=SIMDIR),
		"{SIMDIR}reduced_GTF/GRCh37.85_chr19_22_reduced_exon.gtf".format(SIMDIR=SIMDIR),
		"{SIMDIR}reduced_GTF/GRCh37.85_chr19_22_reduced_me_exon.gtf".format(SIMDIR=SIMDIR)


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
	output: "{SIMDIR}{SAMPLENAME}.isoforms.results"
	threads: CORES
	shell:
		"rsem-calculate-expression --paired-end {input.fa1} {input.fa2} {RSEMREF} {SIMDIR}{SAMPLENAME} --strandedness reverse -p {CORES} --star --star-gzipped-read-file"
		#--fragment-lenth-min 200   --> should we set this parameter?


### modify the RSEM model file (remove the probability of generating reads with quality 2 to prevent sampling of reads with overall low quality)
rule parse_markov_prob:
	input: "{SIMDIR}{SAMPLENAME}.stat/{SAMPLENAME}.model"
	output: "{SIMDIR}{SAMPLENAME}.stat/{SAMPLENAME}_markov_prob"
	shell: "sed -n '14,114p;115q' {input} > {output}"

rule modify_markov_prob:
	input: "{SIMDIR}{SAMPLENAME}.stat/{SAMPLENAME}_markov_prob"
	output: "{SIMDIR}{SAMPLENAME}.stat/{SAMPLENAME}_markov_prob_modified"
	script: "scripts/modify_RSEM_quality_score_distribution.R"

rule modify_RSEM_model:
	input:
		model="{SIMDIR}{SAMPLENAME}.stat/{SAMPLENAME}.model",
		modified_markov_prob = "{SIMDIR}{SAMPLENAME}.stat/{SAMPLENAME}_markov_prob_modified"
	output: "{SIMDIR}{SAMPLENAME}.stat/{SAMPLENAME}.model_modified"
	shell:
		"sed -n '1,13p;14q' {input.model} > {output} \n"
		"cat {input.modified_markov_prob} >> {output} \n"
		"sed -n '115,$p' {input} >> {output}"

rule simulate_data:
	input: model = "{SIMDIR}{SAMPLENAME}.stat/{SAMPLENAME}.model_modified".format(SIMDIR=SIMDIR, SAMPLENAME=SAMPLENAME),
		iso = "{SIMDIR}{SAMPLENAME}.isoforms.results".format(SIMDIR=SIMDIR,SAMPLENAME=SAMPLENAME)
	output: protected( expand("{SIMDIR}simulated_data/simulated_reads_chr19_22_{nr}.fq", SIMDIR=SIMDIR, nr = [1,2]) )
	shell:
		"rsem-simulate-reads {RSEMREF} {input.model} {input.iso} {THETA} {N_READS} {SIMDIR}simulated_data/simulated_reads_chr19_22 --seed {SEED}"


##################

# Count exon truth

##################
rule exon_truth:
	input:
		gtf = GTF,
		sim_iso_res = "{SIMDIR}simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
		fastq1 = "{SIMDIR}simulated_data/simulated_reads_chr19_22_1.fq",
		fastq2 = "{SIMDIR}simulated_data/simulated_reads_chr19_22_2.fq"
	output:
		"{SIMDIR}analysis/GRCh37.85_chr19_22_exon_truth.txt"
	threads: CORES
	script:
		"scripts/count_exon_truth.R"

##################

# Reduce GTF annotation

##################
rule reduce_GTF:
	input:
		gtf=GTF,
		sim_iso_res = "{SIMDIR}simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
	output:
		me = "{SIMDIR}reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf",
		exon = "{SIMDIR}reduced_GTF/GRCh37.85_chr19_22_reduced_exon.gtf",
		me_exon = "{SIMDIR}reduced_GTF/GRCh37.85_chr19_22_reduced_me_exon.gtf"
	script:
		"scripts/reduce_GTF.R"

##################

# Quality control

##################
rule run_fastqc:
	input:
		fastq1 = "{SIMDIR}simulated_data/simulated_reads_chr19_22_1.fq",
		fastq2 = "{SIMDIR}simulated_data/simulated_reads_chr19_22_2.fq"
	output:
		"{SIMDIR}fastqc/simulated_reads_chr19_22_1_fastqc.html",
		"{SIMDIR}fastqc/simulated_reads_chr19_22_2_fastqc.html"
	shell:
		"fastqc -o {SIMDIR}fastqc/ {input.fastq1} {input.fastq2} "
