## Rules to run the R scripts that create the comparison ggplots.

# Count number of reads that overlap with each exon, 

## But it is hard to know what the correct way to do this is. Count all reads taht overlap with an exon? Count only the reads that were actually simulated from this exon? But then we might discard reads that are wrongly mapped to the exon...


# what about only looking at the novel exon? See how many reads were mapped to the novel exons (only count the reads that are actually simulated from these exons) 
# This way we ignore all annotated exons, but we assume that the mappers perform similarly and this is not the focus of our study. 
# We are only interested in the novel exons, because we want to know if the different mappers manage to map reads across the splice junctions, maybe dependent on the actual expression of the novel exon. If the exon is higher expresed, we expect that the mapping will be much better than for a lowly expressed exon.

# --> plot the true novel exon count together with the exon count of the different mappers?

import os


rule compute_mapped_truth_hisat:
    input:
        "Rout/R_packages_install_state.txt",
        bam = "simulation/mapping/{mapper}/{which_reduced_gtf}/{bam_name}_s.bam",
        gtf = config["gtf"],
        sim_iso_res = "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
        removed_gtf = lambda wildcards: config["reduced_exons"][wildcards.which_reduced_gtf]
    params:
        outprefix = "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{bam_name}_"
    conda:
        "../envs/r_scripts.yaml"
    output:
        expand("simulation/mapped_truth/{{mapper}}/{{which_reduced_gtf}}/{{bam_name}}{suffix}", suffix=["_offset_counts.txt", "_offset_counts_removed_exons.txt", "_offset_soft_clipped.txt", "_offset_soft_clipped_removed_exons.txt"])
    script:
        "../scripts/mapped_truth.R"


rule compute_mapped_truth_star:
    input:
        "Rout/R_packages_install_state.txt",
        bam = "simulation/mapping/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam",
        gtf = config["gtf"],
        sim_iso_res = "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
        removed_gtf = lambda wildcards: config["reduced_exons"][wildcards.which_reduced_gtf]
    params:
        outprefix = "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}_"
    conda:
        "../envs/r_scripts.yaml"
    output: 
        expand("simulation/mapped_truth/{{mapper}}/{{which_reduced_gtf}}/{{test_dirnames}}/{{bam_name}}{suffix}", suffix=["_offset_counts.txt", "_offset_counts_removed_exons.txt", "_offset_soft_clipped.txt", "_offset_soft_clipped_removed_exons.txt"])
    script:
        "../scripts/mapped_truth.R"


rule plot_mapping_offsets:
    input:
        "Rout/R_packages_install_state.txt",
        expand("simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{bam_name}{suffix}", mapper = "hisat2", which_reduced_gtf = config["reduced_gtf"], bam_name = "hisat2", suffix = ["_offset_counts.txt", "_offset_counts_removed_exons.txt", "_offset_soft_clipped.txt", "_offset_soft_clipped_removed_exons.txt"]),
        expand("simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}{suffix}", mapper = "STAR", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"], bam_name = "pass2_Aligned.out", suffix = ["_offset_counts.txt", "_offset_counts_removed_exons.txt", "_offset_soft_clipped.txt", "_offset_soft_clipped_removed_exons.txt"]),
         offset_dir = "simulation/mapped_truth"
    params:
        outdir = "simulation/analysis/mapped_offset"
    conda:
        "../envs/r_scripts.yaml"
    output:
        expand("simulation/analysis/mapped_offset/{prefix}{suffix}", prefix = ["all_reads", "reads_removed_exons"], suffix = ["_read_offset_table.txt", "_offset_without_sc.txt"])
    script:
        "../scripts/plot_mapping_offset.R"

###############
## Quality scores at different positions relative to start of soft-clipped regions
rule plot_quality_scores_star:
    input:
        "Rout/R_packages_install_state.txt",
        bam ="simulation/mapping/{mapper}/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_s.bam"
    params:
        outprefix = "simulation/analysis/mapped_offset/sc_quality_score/{mapper}/{which_reduced_gtf}_{test_dirnames}",
        title = "{mapper}: {which_reduced_gtf}; {test_dirnames}"
    output:
        "simulation/analysis/mapped_offset/sc_quality_score/{mapper}/{which_reduced_gtf}_{test_dirnames}_quality_scores_per_position.pdf"
    conda:
        "../envs/r_scripts.yaml"
    script:
        "../scripts/sc_quality.R"

rule plot_quality_scores_hisat:
    input:
        "Rout/R_packages_install_state.txt",
        bam ="simulation/mapping/{mapper}/{which_reduced_gtf}/hisat2_s.bam"
    params:
        outprefix = "simulation/analysis/mapped_offset/sc_quality_score/{mapper}/{which_reduced_gtf}",
        title = "{mapper}: {which_reduced_gtf}"
    output:
        "simulation/analysis/mapped_offset/sc_quality_score/{mapper}/{which_reduced_gtf}_quality_scores_per_position.pdf"
    conda:
        "../envs/r_scripts.yaml"
    script:
        "../scripts/sc_quality.R"


rule plot_quality_scores_star_real_data:
    input:
        "Rout/R_packages_install_state.txt",
        bam ="SRR3192428/mapping/{mapper}/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out_chr19_22_s.bam"
    params:
        outprefix = "SRR3192428/analysis/sc_quality_score/{mapper}/{which_reduced_gtf}_{test_dirnames}",
        title = "{mapper}: {which_reduced_gtf}; {test_dirnames}"
    output:
        "SRR3192428/analysis/sc_quality_score/{mapper}/{which_reduced_gtf}_{test_dirnames}_quality_scores_per_position.pdf"
    conda:
        "../envs/r_scripts.yaml"
    script:
        "../scripts/sc_quality.R"

rule plot_quality_scores_hisat_real_data:
    input:
        "Rout/R_packages_install_state.txt",
        bam ="SRR3192428/mapping/{mapper}/{which_reduced_gtf}/hisat2_s.bam"
    params:
        outprefix = "SRR3192428/analysis/sc_quality_score/{mapper}/{which_reduced_gtf}",
        title = "{mapper}: {which_reduced_gtf}"
    output:
        "SRR3192428/analysis/sc_quality_score/{mapper}/{which_reduced_gtf}_quality_scores_per_position.pdf"
    conda:
        "../envs/r_scripts.yaml"
    script:
        "../scripts/sc_quality.R"

## -----------------------------------------------

rule mapped_truth_sj_hisat:
    input:
        "Rout/R_packages_install_state.txt",
        bam = "simulation/mapping/{mapper}/{which_reduced_gtf}/{bam_name}_s.bam",
        gtf = config["gtf"],
        sim_iso_res = "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
        removed_gtf = lambda wildcards: config["reduced_exons"][wildcards.which_reduced_gtf]
    params:
        outprefix = "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{bam_name}"
    conda:
        "../envs/R_3.5.1.yaml"
    output:
        expand("simulation/mapped_truth/{{mapper}}/{{which_reduced_gtf}}/{{bam_name}}{suffix}", suffix=["_evaluation_SJ_all.txt", "_evaluation_SJ_overl_removed_exons.txt"])
    script:
        "../scripts/mapped_truth_with_sj.R"


rule compute_mapped_truth_sj_star:
    input:
        "Rout/R_packages_install_state.txt",
        bam = "simulation/mapping/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}_s.bam",
        gtf = config["gtf"],
        sim_iso_res = "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results",
        removed_gtf = lambda wildcards: config["reduced_exons"][wildcards.which_reduced_gtf]
    params:
        outprefix = "simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/{bam_name}"
    conda:
        "../envs/R_3.5.1.yaml"
    output: 
        expand("simulation/mapped_truth/{{mapper}}/{{which_reduced_gtf}}/{{test_dirnames}}/{{bam_name}}{suffix}", suffix=["_evaluation_SJ_all.txt", "_evaluation_SJ_overl_removed_exons.txt"])
    script:
        "../scripts/mapped_truth_with_sj.R"


rule plot_mapped_truth_sj:
    input:
        "Rout/R_packages_install_state.txt",
        expand("simulation/mapped_truth/{mapper}/{which_reduced_gtf}/hisat2{suffix}", mapper="hisat2", which_reduced_gtf = config["reduced_gtf"], suffix=["_evaluation_SJ_all.txt", "_evaluation_SJ_overl_removed_exons.txt"]),
        expand("simulation/mapped_truth/{mapper}/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out{suffix}",  mapper = "STAR", which_reduced_gtf = config["reduced_gtf"], test_dirnames = config["star_param"],  suffix=["_evaluation_SJ_all.txt", "_evaluation_SJ_overl_removed_exons.txt"]),
        indir = "simulation/mapped_truth"
    output:
        expand("simulation/analysis/mapped_sj_eval/{prefix}{suffix}", prefix = ["all_reads", "reads_removed_exons"], suffix = ["_evaluation_SJ_barplot.pdf", "_evaluation_SJ_barplot_percent.pdf", "_BAM_unique_mapped_barplot.pdf", "_evaluation_SJ_barplot_accuracy.pdf", "_evaluation_SJ_ROC.pdf", "_evaluation_SJ_PR.pdf", "_evaluation_SJ_barplot_F1.pdf"])
    params:
        outdir = "simulation/analysis/mapped_sj_eval"
    conda:
        "../envs/r_scripts.yaml"
    script:
        "../scripts/plot_mapped_truth_sj.R"

