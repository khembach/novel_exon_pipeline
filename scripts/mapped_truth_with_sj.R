## Compare the mapped and the true location of a read. The mapped coordinates
## including splice junctions are considered.

library(rtracklayer)
library(GenomicAlignments)
library(stringr)
library(dplyr)
library(GenomicFeatures)
library(data.table)

BAM <- snakemake@input[["bam"]]
GTF <- snakemake@input[["gtf"]]
SIM_ISOFORMS_RESULTS <- snakemake@input[["sim_iso_res"]]
OUTPREFIX <- snakemake@params[["outprefix"]]
REMOVED_GTF <- snakemake@input[["removed_gtf"]]


# GTF <- "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
# SIM_ISOFORMS_RESULTS <- "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results"
# BAM <- "simulation/mapping/STAR/me_exon/default/pass2_Aligned.out_s.bam"
# REMOVED_GTF <- "simulation/reduced_GTF/removed_me_exon.gtf"
# OUTPREFIX <- "simulation/mapped_truth/star/me_exon/default/pass2_Aligned.out"


#' Exon exon junction coordinates
#'
#' Compute the transcriptomic coordinates of exon exon junctions
#'
#' @param u integer vector odered widths of all exons in the transcript
#'
#' @return IRanges object with exon-exon junctions
#' @export
#'
#' @examples
get_exon_junction <- function(u){
  cs <- cumsum(u)
  return(IRanges(start=cs[-length(cs)], end=cs[-length(cs)]+1))
}


#' Evaluate gaps in aligned reads
#'
#' Comparison of read gaps with the true set of gaps in the data set. 
#'
#' @param q1 GRanges object with gaps per first reads
#' @param q2 GRanges object with gaps per second reads 
#' @param t1 GRanges object with eej per first reads 
#' @param t2 GRanges object with eej per second reads
#' @param aln GAlignmentPairs with all mapped read
#'
#' @return data.frame with TP, FP, TN, FN of first and second reads
#' @export
#'
#' @examples
evaluat_read_sj <- function(q1, q2, t1, t2, aln){
   ## Identify wrong gaps by comparing the rows between the two data.frames

  ## FP: all reads that do not have a gap in the truth or a wrong gap
  fp1 <- dplyr::setdiff(q1, t1 ) %>% 
    dplyr::pull(read_id) %>% 
    unique %>% length
  fp2 <- dplyr::setdiff(q2, t2 ) %>% 
    dplyr::pull(read_id) %>%
    unique %>% length
  ## TP: all reads with correct gaps
  tp1 <- length(unique(q1$read_id)) - fp1
  tp2 <- length(unique(q2$read_id)) - fp2
  ## TN: all reads that do not have a gap in both the mapping and the truth
  r1_no_gap <- names(aln)[!names(aln) %in% q1$read_id]
  r2_no_gap <- names(aln)[!names(aln) %in% q2$read_id]
  tn1 <- length(r1_no_gap[!r1_no_gap %in% t1$read_id])
  tn2 <- length(r2_no_gap[!r2_no_gap %in% t2$read_id])
  ## FN: all reads that are mapped without a gap but have a gap in the truth
  fn1 <- sum(r1_no_gap %in% t1$read_id)
  fn2 <- sum(r2_no_gap %in% t2$read_id)
  
  data.frame(measure = c("TP", "TP", "FP", "FP", "TN", "TN", "FN", "FN"),
             count = c(tp1, tp2, fp1, fp2, tn1, tn2, fn1, fn2),
             read = rep(c("first", "second"), 4))
}



##------------------------------------------------------------------------------


print("Constructing transcript ranges of all reads")

gtf <- import(GTF)

iso_results <- read.table(SIM_ISOFORMS_RESULTS, header=TRUE)
iso_results <- cbind(iso_results,sid=1:nrow(iso_results))
## sid = id that represents which transcript this read is simulated from

aln <- readGAlignmentPairs(BAM, strandMode=2, use.names=TRUE)
## get the transcriptomic ranges from each simulated read
read_tr_range <-  str_split(string = names(aln), pattern = "_", 
                            simplify=TRUE)[,c(3, 4, 5)] 
colnames(read_tr_range) <- c("sid", "pos", "insertL")
read_tr_range <- apply(read_tr_range,  2, as.numeric)
read_tr_range <- as.data.frame(read_tr_range)
read_tr_range$read_id <- names(aln)

read_tr_range <- read_tr_range %>% 
  dplyr::left_join(dplyr::select(iso_results, transcript_id, length, sid), 
                   by="sid")

read_tr_range <- read_tr_range %>% 
  dplyr::mutate(start1 = length-pos-100,  ## read length 101
                end1 = length-pos,
                start2 = length-pos-insertL+1,
                end2 = length-pos-insertL+101
  )

tr_ranges1 <- data.frame(start = read_tr_range$start1, 
                         end = read_tr_range$end1,
                         transcript_id = as.character(read_tr_range$transcript_id),
                         read_id = read_tr_range$read_id,
                         stringsAsFactors = FALSE)

tr_ranges2 <- data.frame(start = read_tr_range$start2, 
                         end = read_tr_range$end2,
                         transcript_id = as.character(read_tr_range$transcript_id),
                         read_id = read_tr_range$read_id,
                         stringsAsFactors = FALSE)


##------------------------------------------------------------------------------
## Map the transcriptomic coordinates to genomic coordinates -------------------
print("Mapping transcriptomic to genomic coordinates")

## we identify all reads that overlap exon-exon boundaries on the transcript
## to do that, we identify the exon-exon boundaries in transcriptomic coordinates
## we take the exon-exon boundary position and map it to genomic coordinates
## we extract the SJ from the "N" positions in the CIGAR string of each read
## we compare the true SJ with the mapped SJ 
## this way, we do not evalues the actual mapped bases, but only the SJ
## so it does not matter if the reads are soft-clipped or not

## List of exons per transcript
txdb <- makeTxDbFromGRanges(gtf)
tr_granges <- exonsBy(txdb, by="tx", use.names=TRUE)
tr <- transcripts(txdb)

## keep all transcripts with more than one exon
tr_granges_sj <- tr_granges[lengths(tr_granges)>1]

## The exons are already ordered correctly (descending for the "-" strand: 5'
## exon is the last) --> the cumsumm of the exon lengths gives us the correct 
## location of the exon exon junction
tr_eej <- lapply(width(tr_granges_sj), get_exon_junction)
tr_eej <- as(tr_eej, "IRangesList")

# ## get the eej per read, together with the read_id
tr_eej <- tr_eej %>% 
  unlist %>%
  as.data.frame %>%
  dplyr::rename(transcript_id = names)

## join the eej with the reads
## only keep the eej that are within the exon boundaries
r1_truth_eej  <- tr_ranges1 %>% 
  dplyr::inner_join(tr_eej, by = "transcript_id") %>%
  dplyr::filter(start.x <= start.y & end.x >= end.y) %>% 
  dplyr::select(start.y, end.y, read_id, transcript_id) %>%
  dplyr::rename(start = start.y, end = end.y, seqnames = transcript_id)
r2_truth_eej  <- tr_ranges2 %>% 
  dplyr::inner_join(tr_eej, by = "transcript_id") %>%
  dplyr::filter(start.x <= start.y & end.x >= end.y) %>% 
  dplyr::select(start.y, end.y, read_id, transcript_id) %>%
  dplyr::rename(start = start.y, end = end.y, seqnames = transcript_id)

## add the strand information
r1_truth_eej$strand <- strand(tr)[match(r1_truth_eej$seqnames, 
                                        mcols(tr)$tx_name)] %>% as.character
  
r2_truth_eej$strand <- strand(tr)[match(r2_truth_eej$seqnames, 
                                        mcols(tr)$tx_name)] %>% as.character


## genomic coordinates of the exon exon junctions (including last and first base of exon)
##  XXXXX---XXXXXX  annotation
##      xxxxx       genomic eej coordinates
read_id1 <- r1_truth_eej$read_id
r1_truth_eej <- mapFromTranscripts(GRanges(r1_truth_eej), tr_granges)
read_id2 <- r2_truth_eej$read_id
r2_truth_eej <- mapFromTranscripts(GRanges(r2_truth_eej), tr_granges)
## add the read name
mcols(r1_truth_eej)$read_id <- read_id1[mcols(r1_truth_eej)$xHits]
mcols(r2_truth_eej)$read_id <- read_id2[mcols(r2_truth_eej)$xHits]

## remove the last and first base of the touching exons from the range, we only 
## want the sj
r1_truth_eej <- narrow(r1_truth_eej, start=2, end=-2)
r2_truth_eej <- narrow(r2_truth_eej, start=2, end=-2)

## convert back to data.frame
r1_truth_eej <- r1_truth_eej %>% 
  data.frame %>%  
  dplyr::select(-c(xHits, transcriptsHits)) %>%
  dplyr::mutate(strand = droplevels(strand))
r2_truth_eej <- r2_truth_eej %>% 
  data.frame %>%  
  dplyr::select(-c(xHits, transcriptsHits)) %>%
  dplyr::mutate(strand = droplevels(strand))


## Filter out all reads without splice Junctions -------------------------------
print("Preparing GRanges of all gaps")
s1 <- GenomicAlignments::first(aln)[grepl("N", cigar(GenomicAlignments::first(aln)))]
s2 <- GenomicAlignments::second(aln)[grepl("N", cigar(GenomicAlignments::second(aln)))]

## Find the location of the "N" in the read: https://support.bioconductor.org/p/75307/
## and add the start location of the read --> genomic location of gap
s1_range <- cigarRangesAlongReferenceSpace(cigar(s1), ops = "N", pos = start(s1))
s2_range <- cigarRangesAlongReferenceSpace(cigar(s2), ops = "N", pos = start(s2))
names(s1_range) <- names(s1)
names(s2_range) <- names(s2)
s1_range <- unlist(s1_range)
s2_range <- unlist(s2_range)

## The strand of the pair is the strand of its last alignment --> revert first
s1_range <- s1_range %>%
  data.frame %>%
  dplyr::rename(read_id = names) %>%
  dplyr::mutate(seqnames = factor(seqnames(s1[match(read_id, names(s1))])),
                strand = factor(strand(
                  GenomicAlignments::second(aln)[match(read_id, names(aln))])))

s2_range <- s2_range %>%
  data.frame %>%
  dplyr::rename(read_id = names) %>%
  dplyr::mutate(seqnames = factor(seqnames(s2[match(read_id, names(s2))])),
                strand = factor(strand(s2[match(read_id, names(s2))])))


##----------------------------------------------------------------------------
## compare the read gaps to the true eej per read
print("Evaluating all reads")

res <- evaluat_read_sj(s1_range, s2_range, r1_truth_eej, r2_truth_eej, aln)
write.table(res, file = file.path(paste0(OUTPREFIX, "_evaluation_SJ_all.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE)



##----------------------------------------------------------------------------
## We want all reads that were simulated from one of the exons that were
## removed from the gtf annotation.
r_gtf <- import(REMOVED_GTF)

## only keep the read pairs where any of the single reads overlaps with the 
## location of the removed exons
print("Only reads from the removed exons")

read_id1 <- tr_ranges1$read_id
read_id2 <- tr_ranges2$read_id
r1_truth <- mapFromTranscripts(
  GRanges(seqnames = tr_ranges1$transcript_id, 
          ranges = IRanges(start = tr_ranges1$start, end = tr_ranges1$end)), 
  tr_granges)
r2_truth <- mapFromTranscripts(
  GRanges(seqnames = tr_ranges2$transcript_id,
          ranges = IRanges(start = tr_ranges2$start, end = tr_ranges2$end)),
  tr_granges)

mcols(r1_truth)$read_id <- read_id1[mcols(r1_truth)$xHits]
mcols(r2_truth)$read_id <- read_id2[mcols(r2_truth)$xHits]
r1_removed <- subsetByOverlaps(r1_truth, r_gtf)
r2_removed <- subsetByOverlaps(r2_truth, r_gtf)

r_removed <- union(mcols(r1_removed)$read_id, mcols(r2_removed)$read_id)

## Filter out all reads without splice Junctions -------------------------------
print("Removed: Preparing GRanges of all gaps")
aln_removed <- aln[r_removed]
s1 <- GenomicAlignments::first(aln_removed)[
  grepl("N", cigar(GenomicAlignments::first(aln_removed)))]
s2 <- GenomicAlignments::second(aln_removed)[
  grepl("N", cigar(GenomicAlignments::second(aln_removed)))]

## Find the location of the "N" in the read: https://support.bioconductor.org/p/75307/
## and add the start location of the read --> genomic location of gap
s1_range <- cigarRangesAlongReferenceSpace(cigar(s1), ops = "N", pos = start(s1))
s2_range <- cigarRangesAlongReferenceSpace(cigar(s2), ops = "N", pos = start(s2))
names(s1_range) <- names(s1)
names(s2_range) <- names(s2)
s1_range <- unlist(s1_range)
s2_range <- unlist(s2_range)

### data.frame
s1_range <- s1_range %>%
  data.frame %>%
  dplyr::rename(read_id = names) %>%
  dplyr::mutate(seqnames = factor(seqnames(s1[match(read_id, names(s1))])),
                strand = factor(strand(
                  GenomicAlignments::second(aln_removed)[match(read_id,
                                                       names(aln_removed))])))
s2_range <- s2_range %>%
  data.frame %>%
  dplyr::rename(read_id = names) %>%
  dplyr::mutate(seqnames = factor(seqnames(s2[match(read_id, names(s2))])),
                strand = factor(strand(s2[match(read_id, names(s2))])))


##----------------------------------------------------------------------------
## compare the read gaps to the true eej per read
print("Removed: Evaluating all reads")

res <- evaluat_read_sj(s1_range, s2_range, r1_truth_eej, r2_truth_eej, 
                       aln_removed)
write.table(res, file = file.path(paste0(OUTPREFIX, "_evaluation_SJ_overl_removed_exons.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE)

##----------------------------------------------------------------------------



## conda env for R-3.5.1 is  /home/Shared/kathi/microexon_pipeline/.snakemake/conda/a0a0dd0e
# R_LIBS=/home/Shared/Rlib/release-3.5-lib/ /usr/local/R/R-3.4.0/bin/R
