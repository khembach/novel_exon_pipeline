#### These rules map the original FASTQ files because we want to compare the qualtiy scores of soft-clipped bases in the simulation and real data

### STAR is really slow, if it maps against a partial genome, because it tries to find alignments for all reads. It is much faster, if we map against the complete genome and extract only the reads that we are interested in. 
## --> we use the full genome fasta and the GTF file for all chromosomes except 19 and 22. We add our GTF file for the two chromosomes.


rule make_GTF_real_fastq:
    input:
        original_gtf = "/home/Shared/data/seq/microexon_simulation/microexon_study/annotation/Homo_sapiens.GRCh37.85.gtf.gz",
        reduced_gtf = get_gtf
    output:
        "SRR3192428/annotation/Homo_sapiens.GRCh37.85_reduced_{which_reduced_gtf}_chr19_22.gtf"
    shell:
        "zcat ""/home/Shared/data/seq/microexon_simulation/microexon_study/annotation/Homo_sapiens.GRCh37.85.gtf.gz"" | grep -vP '^19\\t|^22\\t' > {output}; "
        "cat {input.reduced_gtf} >> {output}"


rule generate_star_index_real_fastq:
    input:
        fasta = config["genome"],
        gtf = "SRR3192428/annotation/Homo_sapiens.GRCh37.85_reduced_{which_reduced_gtf}_chr19_22.gtf"
    output:
        "SRR3192428/reference/STAR/{which_reduced_gtf}/Genome",
        outdir = "SRR3192428/reference/STAR/{which_reduced_gtf}/"
    threads: 
        config["cores"]
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.outdir} --genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} -sjdbOverhang 100"


rule star_mapping_real_fastq:
    input:
        fastq1 = config["FASTQ1"],
        fastq2 = config["FASTQ2"],
        star_index = "SRR3192428/reference/STAR/{which_reduced_gtf}/"
    output:
        bam1 = temp("SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/Aligned.out.bam"),
        bam2 = temp("SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_Aligned.out.bam"),
        sj1 = "SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/SJ.out.tab",
        sj2 = "SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}/pass2_SJ.out.tab"
    params:
        outdir = "SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}",
        test_param = get_star_param
    threads:
        config["cores"]
    log:
        log1 = "logs/SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}.log",
        log2 = "logs/SRR3192428/mapping/STAR/{which_reduced_gtf}/{test_dirnames}_pass2.log"
    shell:
        "STAR --genomeDir {input.star_index} --readFilesIn {input.fastq1}"
        " {input.fastq2} --readFilesCommand zcat --runThreadN {threads}"
        " --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted"
        " --outFileNamePrefix {params.outdir}/ --sjdbOverhang 100"
        " --outSAMstrandField intronMotif --outSJfilterDistToOtherSJmin"
        " 10 0 0 0 {params.test_param} 2> {log.log1}; "
        "STAR --genomeDir {input.star_index} --readFilesIn {input.fastq1}"
        " {input.fastq2} --readFilesCommand zcat --runThreadN {threads}"
        " --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted"
        " --outFileNamePrefix {params.outdir}/pass2_ --sjdbOverhang 100"
        " --outSAMstrandField intronMotif --outSJfilterDistToOtherSJmin"
        " 10 0 0 0 --sjdbFileChrStartEnd {output.sj1} {params.test_param}"
        " 2> {log.log2}"


rule extract_chr19_22_bam:
    input:
        "SRR3192428/mapping/STAR/{which_reduced_gtf}/{bam_name}_s.bam.bai",
        bam = "SRR3192428/mapping/STAR/{which_reduced_gtf}/{bam_name}_s.bam"
    output:
        "SRR3192428/mapping/STAR/{which_reduced_gtf}/{bam_name}_chr19_22_s.bam"
    shell:
        "samtools view -b {input.bam} 19 22 > {output}"


rule hisat2_mapping_real_fastq:
    input:
       "reference/hisat2/chr19_22/{which_reduced_gtf}/{which_reduced_gtf}_GRCh37.85_chr19_22.1.ht2",
        fastq1 = config["FASTQ1"],
        fastq2 = config["FASTQ2"]
    output:
        sam = temp("SRR3192428/mapping/hisat2/{which_reduced_gtf}/hisat2.sam"),
        novel_splicesites = "SRR3192428/mapping/hisat2/{which_reduced_gtf}/novel_splice_sites_hisat2.txt"
    params:
        basename = "reference/hisat2/chr19_22/{which_reduced_gtf}/{which_reduced_gtf}_GRCh37.85_chr19_22",
        seed = config["SEED"]
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
