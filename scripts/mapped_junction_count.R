## This script computes the number of reads over all junctions that connect to an exon.

library(GenomicAlignments)


REMOVED <- snakemake@input[["removed"]]
BAM <- snakemake@input[["bam"]]
OUTFILE <- snakemake@output[["outfile"]]

## example files
# REMOVED <- "simulation/reduced_GTF/removed_me_exon_unique_classified.txt"
# BAM <- "simulation/mapping/STAR/me_exon/outSJfilterOverhangMin6/Aligned.out_s.bam"


removed <- read.table(REMOVED, header=TRUE) 
removed$count_lsj <- 0  ## if the junction does not exist (terminal exon), the value is NA, otherwise it is the number of reads crossing the junction
removed$count_rsj <- 0
 
# reads <- readGAlignments( BAM, index = BAM, with.which_label=TRUE, param = ScanBamParam(which=GRanges(removed), what = c("qname")) )
# reads <- reads[njunc(reads) > 0,]
# junc <- junctions(reads, use.mcols=TRUE) ## list with all junctions from each reads

## our data comes from Illumina HiSeq 2000 (stranded TruSeq libary preparation with dUTPs ) --> strandMode = 2
readsPair <- readGAlignmentPairs( BAM, index = BAM, with.which_label=TRUE, param = ScanBamParam(which=GRanges(removed), what = c("qname")), strandMode = 2  )
junc <- summarizeJunctions(readsPair)  ## unique list of all junctions


left_ind <- which( !is.na(removed$lend) )
removed[-left_ind, "count_lsj"] <- NA  ## there is no left sj (terminal exon)
lsj <-  GRanges(seqnames = removed$seqnames[left_ind], range=IRanges(removed$lend[left_ind]+1, removed$start[left_ind]-1 ), strand = removed$strand[left_ind] ) ## all left junctions

right_ind <- which(!is.na(removed$rstart))
removed[-right_ind, "count_lsj"] <- NA  ## there is no right sj
rsj <-  GRanges(seqnames = removed$seqnames[right_ind], range=IRanges(removed$end[right_ind]+1, removed$rstart[right_ind]-1 ), strand = removed$strand[right_ind] ) ## all right junctions


m <- match(lsj, junc, ignore.strand=TRUE) 
s <- ifelse(removed[ left_ind[!is.na(m)], "strand" ] == "+", 2, 3)  ## 2 = "plus_score", 3 = "minus_score"; the number of the column with counts that we need for each junction
removed[left_ind[!is.na(m)], "count_lsj" ] <- as.data.frame( mcols(junc[m[!is.na(m)]]) ) [cbind(seq_along(s), s)]

m <- match(rsj, junc, ignore.strand=TRUE) 
s <- ifelse(removed[ right_ind[!is.na(m)], "strand" ] == "+", 2, 3)  ## 2 = "plus_score", 3 = "minus_score"; the number of the column with counts that we need for each junction
removed[right_ind[!is.na(m)], "count_rsj" ] <- as.data.frame( mcols(junc[m[!is.na(m)]]) ) [cbind(seq_along(s), s)]


write.table(removed, OUTFILE, sep = "\t", row.names = FALSE, quote=FALSE )

