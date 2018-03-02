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
        "../scripts/scalmon_exon_counts_scaled.R"




###################

# Derived Salmon counts

###################

# rule derived_Salmon_counts:
#     input:
#             gtf = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf",
#             truth = ,
#             quant = "simulation/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}/quant.sf",
#             fldgz =
#
#
#
#     output:
#         "simulation/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}/quant.sf",
#         outdir = "simulation/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}/"
#     log:
#         "logs/quantification/Salmon/{which_reduced_gtf}/{test_dirnames}.log"
#     script:
#         "salmon quant -i {input.index} -l A -1 {input.fastq1} -2 {input.fastq2} -o {output.outdir} --threads {threads} 2> {log}"









## Rules for the exon quantification with EQP and featureCounts
rule EQP_setup:
    input:
        gtf = GTF,
    output:
        config["eqp_setup"]
    shell:
        "{eqp_dir}/eqp-setup.sh {input.gtf} {ouput}"
