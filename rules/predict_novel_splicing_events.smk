#Predict novel splicing events using the STAR SJ.out.tab file and compare them to the true novel exons that were removed from the reduced GTF

### using SGseq: function analyzeFeatures()
#Splice junctions and exons are predicted from BAM files with predictTxFeatures.
# Known features can be provided as TxFeatures or SGFeatures via argument features.
# If features is not NULL and predict is TRUE, known features are augmented with predictions

###################

# predict novel exons using the junctions from STAR

###################

rule filter_novel_SJ:
    input:
        gtf = lambda wildcards: config["reduced_gtf"][wildcards.which_reduced_gtf],
        sj = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_SJ.out.tab",
        bam = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s.bam"
    output:
        # outfile = "simulation/analysis/filtered_SJ/novel_exons_reduced_{which_reduced_gtf}_{test_dirnames}.txt"
        outfile = "simulation/analysis/filtered_SJ/two_junc_reads/{which_reduced_gtf}/novel_exons_{test_dirnames}.txt"
    log:
        "logs/filter_SJ/{which_reduced_gtf}_{test_dirnames}.log",
    script:
        "../scripts/filter_novel_SJ.R"


###################

# add predicted exons to the GTF file

###################

rule add_predicted_exon_to_gtf:
    input:
        # me = config["reduced_gtf"]["me"],
        gtf = get_gtf,
        exon_prediction = "simulation/analysis/filtered_SJ/{which_reduced_gtf}/novel_exons_{test_dirnames}.txt"
    output:
        outfile = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf"
    script:
        "../scripts/add_novel_exon_to_gtf.R"



###################

# create transcriptome fasta file from the new gtf file

###################
## gffread: http://ccb.jhu.edu/software/stringtie/gff.shtml
rule transcriptome_fasta:
    input:
        gtf = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf",
        genome = config["genome"]
    output:
        outfile = "simulation/transcriptome/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.fasta"
    log:
        "logs/transcriptome_fasta/{which_reduced_gtf}/{test_dirnames}.log"
    script:
        "../scripts/gtf_to_transcript_fasta.R"

# ## alterantive solution with gffread
# ## gffread: http://ccb.jhu.edu/software/stringtie/gff.shtml
# rule transcriptome_fasta:
#     input:
#         gtf = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf",
#         genome = config["genome"]
#     output:
#         outfile = "simulation/transcriptome/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.fasta"
#     log:
#         "logs/transcriptome_fasta/{which_reduced_gtf}/{test_dirnames}.log"
#     shell:
#         "gffread {input.gtf} -g {input.genome} -w {output.outfile}"




###################

# Plot precision-recall curves for the predictions

###################
## output contains last plot to make sure that the script run through
rule plot_PR_curve:
    input:
        removed = "simulation/analysis/mapped_junction_count/removed_{which_reduced_gtf}_unique_classified_{test_dirnames}_junc_count.txt",
        prediction = "simulation/analysis/filtered_SJ/two_junc_reads/{which_reduced_gtf}/novel_exons_{test_dirnames}.txt"
    output:
        "simulation/analysis/exon_prediction_performance/PR/two_junc_reads/{which_reduced_gtf}/{test_dirnames}/PR_expression.png",
        outdir = "simulation/analysis/exon_prediction_performance/PR/two_junc_reads/{which_reduced_gtf}/{test_dirnames}/"
    script:
        "../scripts/plot_PR_curve_novel_sj.R"
