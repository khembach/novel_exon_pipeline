###### This script runs the Rsubread featureCounts on a bam file and saves an R object

## run featurecounts on all BAM files with the different parameters --> easier to compare

GTF <- snakemake@input[["gtf"]]
BAM <- snakemake@input[["bam"]]
OUTFILE <- snakemake@output[["outfile"]]


library(Rsubread)

### We want to summarize the reads on the exon level and not on the metafeature level (e.g. gene)
### we allow reads to overlap with more than one feature (exon), because we want to quantify short microexons and thus expect a read to cover more than one exon!
### 
fc <- featureCounts(files=BAM, annot.ext=GTF, isGTFAnnotationFile = TRUE, GTF.featureType = "exon", useMetaFeatures = FALSE, allowMultiOverlap=TRUE, isPairedEnd = TRUE)
saveRDS(fc, file = OUTFILE)
