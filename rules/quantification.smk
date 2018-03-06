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







#
#
# ###################
#
# # EQP setup
#
# ###################
# ## Rules for the exon quantification with EQP and featureCounts
# rule EQP_setup:
#     input:
#         gtf = GTF,
#     output:
#         config["eqp_setup"]
#     shell:
#         expand("{eqp_dir}/eqp-setup.sh {input.gtf} {output}", eqp_dir = config["eqp"])
#
#
#
# ###################
#
# # EQP-QM exon quantification
#
# ###################
#
# rule EQP_exon_quantification:
#     input:
#         EQP_setup = config["eqp_setup"],
#         bam = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s.bam"
#     output:
#     shell:
#         expand("{eqp_dir}/eqp-quantify.sh -E 3 -J 3 --unambig --unweighted -d {input.EQP_setup} {output.outdir} {input.bam}", eqp_dir = config["eqp"])
