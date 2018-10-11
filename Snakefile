#######
## The directory contains a .Renviron file that specifies the R library. This is necessary, because Rscript does not find the libraries without.

## Or, you can speficy the R, Rscript version and the library in the snakemake file by using the bash login shell:
### use bash login shell for execution
### ATTENTION: snakemake does not load conda environments if you do this!!
# shell.executable("/bin/bash")
# shell.prefix("source ~/.bashrc; ")

### force reexecution of everything in rule all that depends on changed files:
# snakemake -n --forcerun $(snakemake --list-input-changes)
# snakemake --forcerun $(snakemake --list-input-changes) --cores 24 --use-conda


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
include: "rules/mapping_comparison.smk"
include: "rules/mapping.smk"
include: "rules/predict_novel_splicing_events.smk"
include: "rules/quantification.smk"



rule all:
    input:
        # expand("simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam.bai", which_reduced_gtf = config["reduced_gtf"],
        # test_dirnames = config["star_param"], bam_name = ["Aligned.out", "pass2_Aligned.out"]),
        # expand("simulation/mapping/tophat/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam.bai", which_reduced_gtf = config["reduced_gtf"],
        # test_dirnames = config["tophat_param"], bam_name = ["accepted_hits"]),
        # expand("simulation/analysis/removed_exon_truth/{removed_exon}_truth.txt", removed_exon = config["reduced_gtf"]),
        # config["reduced_exons"]["me"],
        # config["reduced_exons"]["exon"],
        # config["reduced_exons"]["me_exon"],
        # config["reduced_gtf"]["me"],
        # config["reduced_gtf"]["exon"],
        # config["reduced_gtf"]["me_exon"]
        # expand("simulation/analysis/removed_exon_truth/removed_{removed_exon}_summary_table.txt", removed_exon = config["reduced_exons"]),
        # expand("simulation/analysis/removed_exon_truth/{removed_exon}_truth.txt", removed_exon = config["reduced_exons"])
        # expand("simulation/analysis/filtered_SJ/{which_reduced_gtf}/novel_exons_{test_dirnames}.txt",
        # which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )
        # expand("simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf",
        # which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )
        # expand("simulation/transcriptome/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.fasta",
        # which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )
        # expand("simulation/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}/quant.sf",
        # which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )
        # # expand("simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam.bai", which_reduced_gtf = "me_exon",
        # test_dirnames = "outSJfilterDistToOtherSJmin0_outSJfilterOverhangMin6", bam_name = ["Aligned.out", "pass2_Aligned.out"])
        # expand("simulation/analysis/derived_Salmon_counts/{which_reduced_gtf}/{test_dirnames}/salmon_coverage_count.txt",
        # which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )
        # expand( "{eqp_setup}/{which_reduced_gtf}/{test_dirnames}/", eqp_setup = config["eqp_setup"],
        # which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )
        # expand("simulation/quantification/EQP/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s-exon.cnt",
        # which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])
        # expand("simulation/quantification/featureCounts/{which_reduced_gtf}/{test_dirnames}/featureCounts.rds",
        # which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])
        # expand( "simulation/analysis/mapped_junction_count/removed_{removed_exon}_unique_classified_{test_dirnames}_junc_count.txt",
        # removed_exon = config["reduced_gtf"], test_dirnames = config["star_param"])
        # expand("simulation/analysis/exon_prediction_performance/PR/{which_reduced_gtf}/{test_dirnames}/PR_expression.png", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"])
        # expand("simulation/analysis/exon_prediction_performance/PR/two_junc_reads_gene_pairs_annotated/{which_reduced_gtf}/{test_dirnames}/PR_expression.pdf", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"] )
        # "simulation/analysis/stringtie/me_exon/outSJfilterOverhangMin6_stringtie.gtf"
        # expand("simulation/analysis/stringtie/PR/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/PR_class_expr.pdf", which_reduced_gtf = config["reduced_gtf"], stringtie_param = config["stringtie_param"], test_dirnames = "outSJfilterOverhangMin6" ),
        # expand("simulation/analysis/stringtie/derived_Salmon_counts/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/salmon_coverage_count.txt",
        # which_reduced_gtf = "me_exon", stringtie_param = config["stringtie_param"], test_dirnames = "outSJfilterOverhangMin6" )
        # "simulation/analysis/stringtie/gffcompare/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6/stringtie.refmap"
        # expand("simulation/analysis/stringtie/gffcompare/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/stringtie.annotated.gtf",
        # which_reduced_gtf = "me_exon", stringtie_param = "minReadCoverage1_minIsoformAbundance0.05", test_dirnames = "outSJfilterOverhangMin6")
        # expand("simulation/analysis/gffcompare/{which_reduced_gtf}/{test_dirnames}/prediction.annotated.gtf",
        # which_reduced_gtf = "me_exon", test_dirnames = "outSJfilterOverhangMin6")
        # expand("simulation/mapping/hisat2/{which_reduced_gtf}/hisat2_s.bam.bai", which_reduced_gtf = config["reduced_gtf"] )
        expand("simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{bam_name}_mapped_truth.txt", mapper = "hisat2", which_reduced_gtf = config["reduced_gtf"], bam_name = "hisat2"),
        expand("simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}_mapped_truth.txt", mapper = "STAR", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"], bam_name = "pass2_Aligned.out")


# rule test:
#     conda:
#         "envs/test.yaml"
#     shell:
#         "python --version"


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
        "simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt"
    threads: CORES
    script:
        "scripts/count_exon_truth.R"

##################

# Reduce GTF annotation

##################
rule reduce_GTF:
    input:
        gtf=GTF,
        truth = "simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt"
    output:
        expand("simulation/reduce_GTF/removed_{removed_exon}_unique.txt", removed_exon = config["reduced_exons"]),
        list(config["reduced_exons"].values()),
        me = config["reduced_gtf"]["me"],
        exon = config["reduced_gtf"]["exon"],
        me_exon = config["reduced_gtf"]["me_exon"]
    script:
        "scripts/reduce_GTF_expressed_exons.R"

rule removed_exons_truth:
    input:
        gtf = lambda wildcards: config["reduced_exons"][wildcards.removed_exon],
        truth = "simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt"
    output:
        outfile = "simulation/analysis/removed_exon_truth/{removed_exon}_truth.txt"
    script:
        "scripts/novel_exon_truth_table.R"

# ## write summary table with information about all removed me/exons
# rule removed_exons_table:
#     input:
#         removed_gtf = lambda wildcards: config["reduced_exons"][wildcards.removed_exon],
#         truth = "simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt",
#         reduced_gtf = lambda wildcards: config["reduced_gtf"][wildcards.removed_exon]
#     output:
#         "simulation/analysis/removed_exon_truth/removed_{removed_exon}_summary_table.txt"
#     script:
#         "scripts/write_removed_exons_table.R"

rule classify_removed_exons:
    input:
        removed = "simulation/reduced_GTF/removed_{removed_exon}_unique.txt"
    output:
        outfile = "simulation/reduced_GTF/removed_{removed_exon}_unique_classified.txt"
    script:
        "scripts/classify_removed_exons.R"


rule mapped_junction_count:
    input:
        removed = "simulation/reduced_GTF/removed_{removed_exon}_unique_classified.txt",
        bam = "simulation/mapping/STAR/{removed_exon}/{test_dirnames}/pass2_Aligned.out_s.bam"
    output:
        outfile = "simulation/analysis/mapped_junction_count/removed_{removed_exon}_unique_classified_{test_dirnames}_junc_count.txt"
    script:
        "scripts/mapped_junction_count.R"

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


##################

# Exon quantification

##################

## featureCounts

## EQP

## Salmon derived counts
## --> discovery of novel splicing events from the STAR splice junctions
## or using SGseq
