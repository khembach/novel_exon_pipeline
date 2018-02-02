## This script takes a gtf file with the novel exons and the truth counts for all exons as input and write a 
## table with only the expression of the novel exons

# GTF <- "simulation/reduced_GTF/removed_me.gtf"
# TRUTH <- "simulation/analysis/GRCh37.85_chr19_22_exon_truth.txt"
# OUTFILE<- "simulation/analysis/removed_exon_truth/me_truth.txt"


GTF <- snakemake@input[["gtf"]]
TRUTH <- snakemake@input[["truth"]]
OUTFILE<- snakemake@output[["outfile"]]


## read the GTF file and truth; overlap the two
## write only the novel exos to file
library(rtracklayer)
gtf <- import(GTF)

truth <- read.table(TRUTH, header=TRUE)
truth_gr <- GRanges(truth)

gtf_truth <- subsetByOverlaps(truth_gr, gtf, type = "equal")

write.table(gtf_truth, file = OUTFILE, quote=FALSE, sep = "\t", row.names = FALSE)
