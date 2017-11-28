### Rules for mapping with STAR and tophat

##################

# Run read mapping

##################
## STAR version STAR_2.5.3a
#### use STAR wrapper!?
### write wrapper for index creation!?


def get_gtf(wildcards):
    return config["reduced_gtf"][wildcards.which_reduced_gtf]  ## e.g. outSJfilterOverhangMin

## generate STAR gemome indices using the different reduced gtf files
rule generate_star_index:
	input:
		# fasta = expand("{GENOMEDIR}{chr}.fa", GENOMEDIR=GENOMEDIR, chr = CHROMS),
		fasta = expand("{GENOMEDIR}Homo_sapiens.GRCh37.dna.chromosome.{chr}.fa", GENOMEDIR=GENOMEDIR, chr = config["chromosomes"]),
		gtf = get_gtf
	output:
		"reference/STAR/chr19_22/{which_reduced_gtf}/Genome",
		outdir = "reference/STAR/chr19_22/{which_reduced_gtf}/"
	threads: CORES
	shell:
		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.outdir} --genomeFastaFiles {input.fasta} "
		"--sjdbGTFfile {input.gtf} -sjdbOverhang 100"


rule generate_bowtie2_index:
	input:
		fasta = expand("{GENOMEDIR}Homo_sapiens.GRCh37.dna.chromosome.{chr}.fa", GENOMEDIR=GENOMEDIR, chr = config["chromosomes"])
		# ",".join("{GENOMEDIR}*.fa")  ## how to join the fasta files with ,?
	output:
		"reference/bowtie2/chr19_22/"
	params:
		prefix = "reference/bowtie2/chr19_22/GRCh37.85_chr19_22"
	threads: CORES
	run:
		fasta_list = ",".join(input.fasta)
		shell("bowtie2-build --seed {SEED} -f {fasta_list} {params.prefix}" )
### use Hisat2 because tophat2 is on low maintenance???


rule mv_fasta_2_bowtie_index:
	input:
		fasta = expand("{GENOMEDIR}Homo_sapiens.GRCh37.dna.chromosome.{chr}.fa", GENOMEDIR=GENOMEDIR, chr = config["chromosomes"]),
		bowtie_index = "reference/bowtie2/chr19_22/"
	output:
		"reference/bowtie2/chr19_22/GRCh37.85_chr19_22.fa"
	shell:
		"cat {input.fasta} >> {output}"


def get_star_param(wildcards):
    return config["star_param"][wildcards.test_dirnames]  ## e.g. outSJfilterOverhangMin
	# return "{SIMDIR}mapping/{test_dirnames}/{{sample}}_pass2.Aligned.out.bam"

### we run 2 pass mapping for all parameters
rule star_mapping:
	input:
		fastq1 = "simulation/simulated_data/simulated_reads_chr19_22_1.fq",
		fastq2 = "simulation/simulated_data/simulated_reads_chr19_22_2.fq",
		star_index = "reference/STAR/chr19_22/{which_reduced_gtf}/"
	output:
		bam1 = temp("simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/Aligned.out.bam"),
		bam2 = temp("simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out.bam"),
		sj1 = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/SJ.out.tab",
		sj2 = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_SJ.out.tab"
		# "{SIMDIR}mapping/{{sample}}.bam"
	params:
		# star_index = "reference/STAR/chr19_22/{which_reduced_gtf}",
		outdir = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}",
		test_param = get_star_param
	threads:
		CORES
	log:
		log1 = "logs/mapping/STAR/{which_reduced_gtf}/{test_dirnames}.log",
		log2 = "logs/mapping/STAR/{which_reduced_gtf}/{test_dirnames}_pass2.log"
	shell:
		"""
		STAR --genomeDir {input.star_index} --readFilesIn {input.fastq1} {input.fastq2} --runThreadN {threads} \
		--outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outFileNamePrefix {params.outdir}/ \
		--sjdbOverhang 100 --outSAMstrandField intronMotif {params.test_param} 2> {log.log1}
		STAR --genomeDir {input.star_index} --readFilesIn {input.fastq1} {input.fastq2} --runThreadN {threads} \
		--outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outFileNamePrefix {params.outdir}/pass2_ \
		--sjdbOverhang 100 --outSAMstrandField intronMotif --sjdbFileChrStartEnd {output.sj1} {params.test_param} 2> {log.log2}
		"""
        ## --outSAMstrandField intronMotif should create an XS attribute in the output file
        ##–outSAMattributes NH HI AS nM XS --> did not work, still no XS attribute in output file!




def get_tophat_param(wildcards):
    return config["tophat_param"][wildcards.test_dirnames]  ## e.g. outSJfilterOverhangMin

rule tophat_mapping:
	input:
		"reference/bowtie2/chr19_22/",
		"reference/bowtie2/chr19_22/GRCh37.85_chr19_22.fa",
		fastq1 = "simulation/simulated_data/simulated_reads_chr19_22_1.fq",
		fastq2 = "simulation/simulated_data/simulated_reads_chr19_22_2.fq",
		gtf = get_gtf
	output:
		bam = temp("simulation/mapping/tophat/{which_reduced_gtf}/{test_dirnames}/accepted_hits.bam")
	params:
		outdir = "simulation/mapping/tophat/{which_reduced_gtf}/{test_dirnames}/",
		test_param = get_tophat_param,
		bowtie2_index = "reference/bowtie2/chr19_22/GRCh37.85_chr19_22"
	threads:
		CORES
	log:
		"logs/mapping/tophat/{which_reduced_gtf}/{test_dirnames}.log"
	conda:
		"../envs/tophat.yaml"  ## use latest tophat and python 2.7
	shell:
		"""
		tophat --num-threads {threads} --max-multihits 1 --GTF {input.gtf} {params.test_param} --output-dir {params.outdir} \
		{params.bowtie2_index} {input.fastq1} {input.fastq2} 2> {log}
		"""


rule sort_bam:
	input:
		bam = "simulation/mapping/{test_dirnames}/{bam_name}.bam",
	output:
		bam_sorted = protected("simulation/mapping/{test_dirnames}/{bam_name}_s.bam"),
	threads:
		4
	shell:
		"samtools sort {input.bam} -o {output.bam_sorted} -@ {threads}"


rule index_bam:
	input:
		bam = "simulation/mapping/{test_dirnames}/{bam_name}_s.bam",
	output:
		bam_indexed = "simulation/mapping/{test_dirnames}/{bam_name}_s.bam.bai",
	threads:
		4
	shell:
		"samtools index -@ {threads} {input.bam}"
