#### These rules map the original FASTQ files because we want to compare the qualtiy scores of soft-clipped bases in the simulation and real data

rule star_mapping_real_fastq:
    input:
        fastq1 = FASTQ1,
        fastq2 = FASTQ2,
        star_index = "reference/STAR/chr19_22/{which_reduced_gtf}/"
    output:
        bam1 = temp("SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/Aligned.out.bam"),
        bam2 = temp("SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out.bam"),
        sj1 = "SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/SJ.out.tab",
        sj2 = "SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_SJ.out.tab"
        # "{SIMDIR}mapping/{{sample}}.bam"
    params:
        outdir = "SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}",
        test_param = get_star_param
    threads:
        config["cores"]
    log:
        log1 = "logs/SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}.log",
        log2 = "logs/SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}_pass2.log"
    shell:
        """
        STAR --genomeDir {input.star_index} --readFilesIn {input.fastq1} {input.fastq2} --runThreadN {threads} \
        --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outFileNamePrefix {params.outdir}/ \
        --sjdbOverhang 100 --outSAMstrandField intronMotif --outSJfilterDistToOtherSJmin 10 0 0 0 {params.test_param} 2> {log.log1}
        STAR --genomeDir {input.star_index} --readFilesIn {input.fastq1} {input.fastq2} --runThreadN {threads} \
        --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outFileNamePrefix {params.outdir}/pass2_ \
        --sjdbOverhang 100 --outSAMstrandField intronMotif \
        --outSJfilterDistToOtherSJmin 10 0 0 0 \
        --sjdbFileChrStartEnd {output.sj1} {params.test_param} 2> {log.log2}
        """

        ## --outSAMstrandField intronMotif should create an XS attribute in the output file
        ##–outSAMattributes NH HI AS nM XS --> did not work, still no XS attribute in output file!


rule hisat2_mapping_real_fastq:
    input:
       "reference/hisat2/chr19_22/{which_reduced_gtf}/{which_reduced_gtf}_GRCh37.85_chr19_22.1.ht2",
        fastq1 = FASTQ1,
        fastq2 = FASTQ2
    output:
        sam = temp("SRR3192428/mapping/hisat2/{which_reduced_gtf}/hisat2.sam"),
        novel_splicesites = "SRR3192428/mapping/hisat2/{which_reduced_gtf}/novel_splice_sites_hisat2.txt"
    params:
        basename = "reference/hisat2/chr19_22/{which_reduced_gtf}/{which_reduced_gtf}_GRCh37.85_chr19_22",
        seed = SEED
    threads:
        config["cores"]
    log:
        "logs/SRR3192428/mapping/hisat2/hisat_{which_reduced_gtf}.log"
    shell:
        """
        hisat2 -x {params.basename} -1 {input.fastq1} -2 {input.fastq2} -S {output.sam} --novel-splicesite-outfile {output.novel_splicesites} --rna-strandness RF -k 1 --new-summary --no-unal -p {threads} --seed {params.seed} 2> {log}
        """

rule convert_bam_real_fastq:
    input:
        sam = "SRR3192428/mapping/hisat2/{which_reduced_gtf}/hisat2.sam"
    output:
        bam = "SRR3192428/mapping/hisat2/{which_reduced_gtf}/hisat2.bam"
    threads:
        4
    shell:
        "samtools view -b -@ {threads} {input.sam} > {output.bam}"


rule sort_bam_real_fastq:
    input:
        bam = "SRR3192428/mapping/{mapper}/{which_reduced_gtf}/{bam_name}.bam",
    output:
        bam_sorted = protected("SRR3192428/mapping/{mapper}/{which_reduced_gtf}/{bam_name}_s.bam"),
    threads:
        4
    shell:
        "samtools sort {input.bam} -o {output.bam_sorted} -@ {threads}"


rule index_bam_real_fastq:
    input:
        bam = "SRR3192428/mapping/{mapper}/{which_reduced_gtf}/{bam_name}_s.bam",
    output:
        bam_indexed = "SRR3192428/mapping/{mapper}/{which_reduced_gtf}/{bam_name}_s.bam.bai",
    threads:
        4
    shell:
        "samtools index -@ {threads} {input.bam}"
