## data was downloaded from ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX160/SRX1603411/SRR3192428/
## GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2072377
FASTQ1 = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/SRR3192428_1.fastq.gz"
FASTQ2 = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/SRR3192428_2.fastq.gz"
GTF = "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
CORES = 16
GENOMEDIR = "/home/Shared/data/seq/microexon_simulation/microexon_study/genome/reference/"
RSEMREF = "reference/RSEM/hg19_chr19_22"
REALSAMPLE = "/home/Shared/data/seq/microexon_simulation/microexon_study/real_data/"
SIMDIR="simulation/"
SAMPLENAME = "SRR3192428"

rule all:
	input: SIMDIR + SAMPLENAME +".isoforms.results"


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
	output: SIMDIR + SAMPLENAME +".isoforms.results"
	threads: CORES
	shell:
		"rsem-calculate-expression --paired-end {input.fa1} {input.fa2} {RSEMREF} {SIMDIR}{SAMPLENAME} --strandedness reverse -p {CORES} --star --star-gzipped-read-file"
		#--fragment-lenth-min 200   --> should we set this parameter?


### modify model file manually (remove the probability of generating reads with quality 2 to prevent sampling of reads with overall low quality)
