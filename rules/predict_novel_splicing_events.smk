#Predict novel splicing events using the STAR SJ.out.tab file and compare them to the true novel exons that were removed from the reduced GTF

### using SGseq: function analyzeFeatures()
#Splice junctions and exons are predicted from BAM files with predictTxFeatures.
# Known features can be provided as TxFeatures or SGFeatures via argument features.
# If features is not NULL and predict is TRUE, known features are augmented with predictions

###################

# predict novel exons using the junctions from STAR

###################

rule predict_novel_exons:
    input:
        script = "R/predict_novel_exons.R",
        gtf = lambda wildcards: config["reduced_gtf"][wildcards.which_reduced_gtf],
        sj = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_SJ.out.tab",
        bam = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s.bam"
    output:
        "simulation/analysis/filtered_SJ/{exon_pred_dir}/{which_reduced_gtf}/novel_exons_{test_dirnames}.txt"
    params:
        overhang_min = lambda wildcards: config["overhang_min"][wildcards.test_dirnames]
    log:
        "logs/Rout/exon_prediction/{which_reduced_gtf}_{test_dirnames}.log"
    shell:
        '''{Rbin} CMD BATCH --no-restore --no-save "--args GTF='{input.gtf}' SJFILE='{input.sj}' BAM='{input.bam}' OVERHANGMIN='{params.overhang_min}' OUTFILE='{output}'" {input.script} {log}'''  


###################

# add predicted exons to the GTF file

###################

rule add_predicted_exon_to_gtf:
    input:
        script = "R/add_predictions_to_gtf.R",
        gtf = get_gtf,
        exon_prediction = "simulation/analysis/filtered_SJ/{exon_pred_dir}/{which_reduced_gtf}/novel_exons_{test_dirnames}.txt"
    output:
        outfile = "simulation/reduced_GTF_with_predicted_exons/{exon_pred_dir}/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf"
    log:
        "logs/Rout/extend_gtf/{which_reduced_gtf}_{test_dirnames}.log"
    shell:
        '''{Rbin} CMD BATCH --no-restore --no-save "--args GTF='{input.gtf}' EXONPRED='{input.exon_prediction}' OUTFILE='{output}'" {input.script} {log}'''


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
        prediction = "simulation/analysis/filtered_SJ/two_junc_reads_gene_pairs_annotated/{which_reduced_gtf}/novel_exons_{test_dirnames}.txt"
    output:
        "simulation/analysis/exon_prediction_performance/PR/two_junc_reads_gene_pairs_annotated/{which_reduced_gtf}/{test_dirnames}/PR_expression.pdf",
        outdir = "simulation/analysis/exon_prediction_performance/PR/two_junc_reads_gene_pairs_annotated/{which_reduced_gtf}/{test_dirnames}/"
    conda:
         "../envs/r_scripts.yaml"
    script:
        "../scripts/plot_PR_curve_novel_sj.R"




###################

# Predict transcripts with StringTie

###################
## the library was prepared with dUTPs, so it is a stranded library fr-firststrand (--rf)
# -a <int>  Junctions that don't have spliced reads that align across them with at least this amount
            # of bases on both sides are filtered out. Default: 10
# -c <float>    Sets the minimum read coverage allowed for the predicted transcripts.
                # A transcript with a lower coverage than this value is not shown in the output. Default: 2.5
#  --> probably we will have to lower this value, most lowly expressed transcripts will not be recovert with 2.5
# -f <0.0-1.0>  Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript
                # assembled at a given locus. Lower abundance transcripts are often artifacts of incompletely spliced precursors
                # of processed transcripts. Default: 0.1
# -j <float>    There should be at least this many spliced reads that align across a junction (i.e. junction coverage). This number
            # can be fractional, since some reads align in more than one place. A read that aligns in n places will contribute 1/n to
            # the junction coverage. Default: 1
# -t    This parameter disables trimming at the ends of the assembled transcripts. By default StringTie adjusts the predicted transcript's start
            # and/or stop coordinates based on sudden drops in coverage of the assembled transcript.
# default
# -a 3
# -a 6 (min value in STAR alignment)
# -c 1
# -f 0.05
# -f 0.2
# -j 1 (STAR value) --> no change required?
# -t (probably we should not set this parameter, but I still want to know what happens if we use it)
#

def get_stringtie_param(wildcards):
    return config["stringtie_param"][wildcards.stringtie_param]  ## e.g. minJuncOverhang6

rule stringtie_assemble_transcripts:
    input:
        gtf = lambda wildcards: config["reduced_gtf"][wildcards.which_reduced_gtf],
        bam = "simulation/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s.bam"
    output:
        "simulation/analysis/stringtie/predictions/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}_stringtie.gtf"
    params:
        test_param = get_stringtie_param
    threads:
        8
    shell:
        "stringtie {input.bam} -o {output} -p {threads} -G {input.gtf} --rf {params.test_param}"


rule stringtie_pred_exons:
    input:
        gtf = lambda wildcards: config["reduced_gtf"][wildcards.which_reduced_gtf],
        strtie = "simulation/analysis/stringtie/predictions/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}_stringtie.gtf"
    output:
        outfile = "simulation/analysis/stringtie/predictions/{which_reduced_gtf}/{stringtie_param}/novel_exons_{test_dirnames}_stringtie.txt"
    script:
        "../scripts/stringtie_novel_exons.R"


rule plot_PR_curve_stringtie:
    input:
        removed = "simulation/analysis/mapped_junction_count/removed_{which_reduced_gtf}_unique_classified_{test_dirnames}_junc_count.txt",
        prediction = "simulation/analysis/stringtie/predictions/{which_reduced_gtf}/{stringtie_param}/novel_exons_{test_dirnames}_stringtie.txt"
    output:
        "simulation/analysis/stringtie/PR/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/PR_class_expr.pdf",
        outdir = "simulation/analysis/stringtie/PR/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/"
    conda:
         "../envs/r_scripts.yaml"
    script:
        "../scripts/plot_PR_curve_novel_sj.R"


rule stringtie_transcriptome_fasta:
    input:
        gtf = "simulation/analysis/stringtie/predictions/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}_stringtie.gtf",
        genome = config["genome"]
    output:
        outfile = "simulation/analysis/stringtie/transcriptome/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}_stringtie_{stringtie_param}.fasta"
    log:
        "logs/stringtie/transcriptome_fasta/{which_reduced_gtf}/{stringtie_param}_{test_dirnames}.log"
    script:
        "../scripts/gtf_to_transcript_fasta.R"



############

## Compare predicted GTF with reference annotation

###########

rule gffcompare_stringtie:
    input:
        stringtie = "simulation/analysis/stringtie/predictions/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}_stringtie.gtf",
        ref_gtf = GTF
    output:
        outfile = "simulation/analysis/stringtie/gffcompare/{which_reduced_gtf}/{stringtie_param}/{test_dirnames}/stringtie.annotated.gtf",
        refmap = "simulation/analysis/stringtie/predictions/{which_reduced_gtf}/{stringtie_param}/stringtie.{test_dirnames}_stringtie.gtf.refmap",
    shell:
        "gffcompare {input.stringtie} -o simulation/analysis/stringtie/gffcompare/{wildcards.which_reduced_gtf}/{wildcards.stringtie_param}/{wildcards.test_dirnames}/stringtie -r {input.ref_gtf} -e 100 -d 100 -V"


rule gffcompare_prediction:
    input:
        prediction = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf",
        ref_gtf = GTF
    output:
        outfile = "simulation/analysis/gffcompare/{which_reduced_gtf}/{test_dirnames}/prediction.annotated.gtf",
        refmap = "simulation/reduced_GTF_with_predicted_exons/{which_reduced_gtf}/prediction.GRCh37.85_chr19_22_novel_exons_{test_dirnames}.gtf.refmap"
    shell:
        "gffcompare {input.prediction} -o simulation/analysis/gffcompare/{wildcards.which_reduced_gtf}/{wildcards.test_dirnames}/prediction -r {input.ref_gtf} -e 100 -d 100 -V"
