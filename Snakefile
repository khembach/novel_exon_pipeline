#######
## The directory contains a .Renviron file that specifies the R library. This is necessary, because Rscript does not find the libraries without.
## Or, you can speficy the R, Rscript version and the library in the snakemake file by using the bash login shell:

### use bash login shell for execution
shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")

## config file
configfile: "config.yaml"


## data was downloaded from ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX160/SRX1603411/SRR3192428/
## GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2072377
FASTQ1 = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/SRR3192428_1.fastq.gz"
FASTQ2 = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/SRR3192428_2.fastq.gz"
GTF = "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
CORES = 12
GENOMEDIR = "/home/Shared/data/seq/microexon_simulation/microexon_study/genome/reference/"
RSEMREF = "reference/RSEM/hg19_chr19_22"
REALSAMPLE = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/"
SAMPLENAME = "SRR3192428"
THETA = 0.05
N_READS = 1118017
## the number of reads that were mapped to chromosome 19 and 22: 1118017
## 8306564 = number of reads = [348875684]/4*{1500+500}/(21000)        [lines in real fastq data] {genes in chr19 + chr22} / (total genes)
SEED = 584

STAR_PARAMS_DIRNAME = ["default", "default_2_pass", "outSJfilterOverhangMin9", "outSJfilterOverhangMin6", "outSJfilterCountTotalMin3",
 "scoreGenomicLengthLog2scale0", "alignSJoverhangMin3"]

### sub sections of the workflow:
include: "rules/rsem_simulation.smk"
include: "rules/mapping.smk"


rule all:
    input:
        #    expand("simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam.bai", which_reduced_gtf = config["reduced_gtf"],
        #    test_dirnames = config["star_param"], bam_name = ["Aligned.out", "pass2_Aligned.out"])
        expand("simulation/mapping/tophat/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam.bai", which_reduced_gtf = config["reduced_gtf"],
        test_dirnames = config["tophat_param"], bam_name = ["accepted_hits"])


##################

# RSEM

##################


##################

# Count exon truth

##################
rule exon_truth:
	input:
		gtf = GTF,
		sim_iso_res = "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
		fastq1 = "simulation/simulated_data/simulated_reads_chr19_22_1.fq",
		fastq2 = "simulation/simulated_data/simulated_reads_chr19_22_2.fq"
	output:
		"simulation/analysis/GRCh37.85_chr19_22_exon_truth.txt"
	threads: CORES
	script:
		"scripts/count_exon_truth.R"

##################

# Reduce GTF annotation

##################
rule reduce_GTF:
    input:
        gtf=GTF,
        sim_iso_res = "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
    output:
        me = config["reduced_gtf"]["me"],
        exon = config["reduced_gtf"]["exon"],
        me_exon = config["reduced_gtf"]["me_exon"]
    script:
        "scripts/reduce_GTF.R"



##################

# Quality control

##################
rule run_fastqc:
	input:
		fastq1 = "simulation/simulated_data/simulated_reads_chr19_22_1.fq",
		fastq2 = "simulation/simulated_data/simulated_reads_chr19_22_2.fq"
	output:
		"simulation/fastqc/simulated_reads_chr19_22_1_fastqc.html",
		"simulation/fastqc/simulated_reads_chr19_22_2_fastqc.html"
	shell:
		"fastqc -o simlation/fastqc/ {input.fastq1} {input.fastq2} "



##################

# Run read mapping

##################
