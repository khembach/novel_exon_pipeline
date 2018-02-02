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
        outfile = "simulation/analysis/filtered_SJ/novel_exons_reduced_{which_reduced_gtf}_{test_dirnames}.txt"
    log:
        "logs/filter_SJ/{which_reduced_gtf}_{test_dirnames}.log",
    script:
        "../scripts/filter_novel_SJ.R"
