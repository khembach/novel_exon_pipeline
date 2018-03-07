###################

# Salmon index generation

###################

rule Salmon_index:
    input:
        fasta  = "simulation/transcriptome/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.fasta"
    output:
        "reference/Salmon/{which_reduced_gtf}/{test_dirnames}/hash.bin",
        dir = "reference/Salmon/{which_reduced_gtf}/{test_dirnames}/"
    shell:
        "salmon index -t {input.fasta} -i {output.dir} --type quasi -k 31"


###################

# Transcript quantification with Salmon

###################

rule Salmon_quant:
    input:
        index = "reference/Salmon/{which_reduced_gtf}/{test_dirnames}/",
        fastq1 = "simulation/simulated_data/simulated_reads_chr19_22_1.fq",
        fastq2 = "simulation/simulated_data/simulated_reads_chr19_22_2.fq"
    output:
        "simulation/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}/quant.sf",
        outdir = "simulation/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}/"
    threads:
        4
    log:
        "logs/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}.log"
    shell:
        "salmon quant -i {input.index} -l A -1 {input.fastq1} -2 {input.fastq2} -o {output.outdir} --threads {threads} 2> {log}"

        # "salmon quant -i {input.index} -l A -1 {input.fastq1} -2 {input.fastq2} -o {output.outdir} --threads {threads} 2> {log}"



###################

# Derived Salmon counts

###################

rule derived_Salmon_counts:
    input:
        gtf = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf",
        truth = "simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt",
        removed = "simulation/reduced_GTF/removed_{which_reduced_gtf}_unique_classified.txt",
        quant = "simulation/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}/quant.sf",
        fldgz = "simulation/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}/aux_info/fld.gz"
    output:
        "simulation/analysis/derived_Salmon_counts/{which_reduced_gtf}/{test_dirnames}/salmon_coverage_count.txt",
        outdir = "simulation/analysis/derived_Salmon_counts/{which_reduced_gtf}/{test_dirnames}/"
    log:
        "logs/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}.log"
    script:
        "../scripts/salmon_exon_counts_scaled.R"




###################

# EQP setup

###################
## Rules for the exon quantification with EQP and featureCounts
rule EQP_setup:
    input:
        gtf = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf"
    output:
        config["eqp_setup"] + "/{which_reduced_gtf}/{test_dirnames}/map-files/GRCh37.85_chr19_22_novel_exons_{test_dirnames}_exon_junction.map.gz",
        setup_dir = config["eqp_setup"] + "/{which_reduced_gtf}/{test_dirnames}"
    conda:
        "/home/Shared/kathi/microexon_pipeline/envs/test.yaml"  ## use python 2.7
    shell:
        config["eqp"] + "/eqp-setup.sh {input.gtf} {output.setup_dir}"


## the environment was installed with
## conda create --name EQP python=2.7 samtools bedtools r numpy argparse gettext
## and the yml file was exported from the activated environment with
## conda env export > EQP.yml


###################

# EQP-QM exon quantification

###################
# -e : compute exon counts (we do not need gene and junction counts!)
# -E int : Minimal overlap of a read with an exon [5]
# -J INT: Minimal overlap of a read with both exons on a junction [8]
# -s STRING: process reads as strand-specific in direction STRING (forward for orientation fr or backward for orientation rf)
# --unambig: count only reads that can be assigned unambiguously to a single gene or exon when creating the gene or exon counts
# --unweighted: do not use read weights in the generation of counts
# --nosort: the alignment file is already sorted by names; do not sort it

rule EQP_exon_quantification:
    input:
        config["eqp_setup"] + "/{which_reduced_gtf}/{test_dirnames}/map-files/GRCh37.85_chr19_22_novel_exons_{test_dirnames}_exon_junction.map.gz",  ## to make sure that the setup completed successfully
        setup_dir = config["eqp_setup"] + "/{which_reduced_gtf}/{test_dirnames}",
        bam = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s.bam"
    output:
        "simulation/quantification/EQP/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s-exon.cnt",
        outdir = "simulation/quantification/EQP/{which_reduced_gtf}/{test_dirnames}/"
    conda:
        "../envs/EQP.yaml"  ## use python 2.7
    shell:
        config["eqp"] + "/eqp-quantify.sh -e -E 3 -J 3 --unambig --unweighted --nosort -d {input.setup_dir} {output.outdir} {input.bam}"


###################

# FeatureCounts exon quantification

###################

rule featureCounts:
    input:
        gtf = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf",
        bam = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s.bam"
    output:
        outfile = "simulation/quantification/featureCounts/{which_reduced_gtf}/{test_dirnames}/featureCounts.rds"
    script:
        "../scripts/run_featureCounts.R"



###################

# Comparison between featureCounts, EQP and the derived Salmon counts

###################
#
# rule quantification_comparison:
#     input:
#     output:
#     script:
#         "scripts/comparison_ggplot_salmon.R"
