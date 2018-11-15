## Compare the mapped and the true location of a read. The mapped coordinates
## including splice junctions are considered.

library(rtracklayer)
library(GenomicAlignments)
library(stringr)
library(dplyr)
# library(ensembldb)

library(GenomicFeatures)
library(data.table)

BAM <- snakemake@input[["bam"]]
GTF <- snakemake@input[["gtf"]]
SIM_ISOFORMS_RESULTS <- snakemake@input[["sim_iso_res"]]
OUTPREFIX <- snakemake@params[["outprefix"]]
REMOVED_GTF <- snakemake@input[["removed_gtf"]]
# ENSEMBL_DB_SQLITE <- snakemake@params[["ensembl_db_sqlite"]]


# GTF <- "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
# SIM_ISOFORMS_RESULTS <- "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results"
# ## OUTFILE <- "simulation/mapped_truth/STAR/me_exon/default_mapped_truth.txt"
# BAM <- "simulation/mapping/STAR/me_exon/default/pass2_Aligned.out_s.bam"
# REMOVED_GTF <- "simulation/reduced_GTF/removed_me_exon.gtf"
# OUTPREFIX <- "simulation/mapped_truth/star/me_exon/default/pass2_Aligned.out"
# ## ENSEMBL_DB_PATH <- "../annotation"
# ## ENSEMBL_DB_NAME <- "Homo_sapiens.GRCh37.85_chr19_22.sqlite"
# ## ENSEMBL_DB_SQLITE <- "annotation/Homo_sapiens.GRCh37.85_chr19_22.sqlite"



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


#' Exon exon junction 
#' 
#' This function creates a data.frame with the seqnames, start, end, strand and
#' read name of all exon exon junctions that overlap the reads from a given
#' transcript.
#'
#' @param tr_id character transcript if
#' @param hits hits object with mappings from read to exon exon junction 
#' @param tr_r_list IRangesList with read locations per transcript
#'
#' @return data.frame that can be converted to GRanges
#' @export
#'
#' @examples
get_eej_df <- function(tr_id, hits, tr_r_list){
  data.frame(seqnames = rep(tr_id, length(hits)), 
             start = start(tr_eej[[tr_id]][subjectHits(hits)]),
             end = end(tr_eej[[tr_id]][subjectHits(hits)]), 
             strand=rep(strand(tr[mcols(tr)$tx_name == tr_id]), length(hits)),
             read_id=mcols(tr_r_list[[tr_id]][queryHits(hits)])$read_id)
}



#' Title
#'
#' @param s1_range_list GRangesList with gaps per reads 1
#' @param s2_range_list GRAngesList with gaps per reads 2
#' @param r1_truth_eej_list GRangesList with eej per reads 1
#' @param r2_truth_eej_list GRangesList with eej per reads 2
#' @param aln GAlignmentPairs with all mapped read
#'
#' @return data.frame with TP, FP, TN, FN of read 1 and read 2
#' @export
#'
#' @examples
evaluat_read_sj_slow <- function(s1_range_list, s2_range_list, 
                            r1_truth_eej_list,r2_truth_eej_list, aln){
  
  ## we only count the read gaps that cannot be matched with a true sj
  ## FP: all reads that do not have a gap in the truth
  fp1 <- which(!names(s1_range_list) %in% names(r1_truth_eej_list))
  fp2 <- which(!names(s2_range_list) %in% names(r2_truth_eej_list))
  
  ## TODO: this is the correct way to evaluate the gaps in all reads!
  ## How do we speed this up?
  ## TP: all reads 
  comp1 <- s1_range_list[-fp1] %in% r1_truth_eej_list[names(s1_range_list[-fp1])]
  comp2 <- s2_range_list[-fp2] %in% r2_truth_eej_list[names(s2_range_list[-fp2])]
  tp1 <- length(comp1[ all(comp1)])
  tp2 <- length(comp2[ all(comp2)])
  ## TN: all reads that do not have a gap in both the mapping and the truth
  r1_no_gap <- names(aln)[!names(aln) %in% names(s1_range_list)]
  r2_no_gap <- names(aln)[!names(aln) %in% names(s2_range_list)]
  tn1 <- length(r1_no_gap[!r1_no_gap %in% names(r1_truth_eej_list)])
  tn2 <- length(r2_no_gap[!r2_no_gap %in% names(r2_truth_eej_list)])
  ## FN: all reads that are mapped without gap but have a gap in the truth
  fn1 <- sum(r1_no_gap %in% names(r1_truth_eej_list))
  fn2 <- sum(r2_no_gap %in% names(r2_truth_eej_list))
  
  data.frame(measure = c("TP", "TP", "FP", "FP", "TN", "TN", "FN", "FN"),
             count = c(tp1, tp2, length(fp1), length(fp2), tn1, tn2, fn1, fn2),
             read = rep(c("first", "second"), 4))
  
}

#' Evaluate gaps in aligned reads
#'
#' Comparison of read gaps with the true set of gaps in the data set. 
#'
#' @param s1_range Granges of 
#' @param s2_range 
#' @param r1_truth_eej 
#' @param r2_truth_eej 
#' @param aln 
#'
#' @return
#' @export
#'
#' @examples
evaluat_read_sj <- function(s1_range, s2_range, 
                            r1_truth_eej, r2_truth_eej, aln){
  ## Convert the GRanges of the read gaps and the true read SJ to data.frame 
  ## identical columns.
  ## Identify wrong gaps by comparing the rows between the two data.frames
  q1 <- s1_range %>% 
    as.data.frame(row.names=NULL, stringsAsFactors=FALSE) %>% 
    dplyr::mutate(read_id = names(s1_range))
  q2 <- s2_range %>% 
    as.data.frame(row.names=NULL, stringsAsFactors=FALSE) %>% 
    dplyr::mutate(read_id = names(s2_range))

  t1 <- r1_truth_eej %>% 
    as.data.frame(stringsAsFactors=FALSE) %>%  
    dplyr::select(-c(xHits, transcriptsHits)) %>%
    dplyr::mutate(read_id = as.character(read_id)) 
  t2 <- r2_truth_eej %>% 
    as.data.frame(stringsAsFactors=FALSE) %>%  
    dplyr::select(-c(xHits, transcriptsHits)) %>%
    dplyr::mutate(read_id = as.character(read_id)) 
  
  ## FP: all reads that do not have a gap in the truth or a wrong gap
  fp1 <- dplyr::setdiff(q1, t1 ) %>% 
    dplyr::pull(read_id) %>% 
    unique %>% length
  fp2 <- dplyr::setdiff(q2, t2 ) %>% 
    dplyr::pull(read_id) %>%
    unique %>% length
  ## TP: all reads with correct gaps
  tp1 <- length(unique(names(s1_range))) - fp1
  tp2 <- length(unique(names(s2_range))) - fp2
  ## TN: all reads that do not have a gap in both the mapping and the truth
  r1_no_gap <- names(aln)[!names(aln) %in% names(s1_range)]
  r2_no_gap <- names(aln)[!names(aln) %in% names(s2_range)]
  tn1 <- length(r1_no_gap[!r1_no_gap %in% mcols(r1_truth_eej)$read_id])
  tn2 <- length(r2_no_gap[!r2_no_gap %in% mcols(r2_truth_eej)$read_id])
  ## FN: all reads that are mapped without a gap but have a gap in the truth
  fn1 <- sum(r1_no_gap %in% mcols(r1_truth_eej)$read_id)
  fn2 <- sum(r2_no_gap %in% mcols(r2_truth_eej)$read_id)
  
  data.frame(measure = c("TP", "TP", "FP", "FP", "TN", "TN", "FN", "FN"),
             count = c(tp1, tp2, fp1, fp2, tn1, tn2, fn1, fn2),
             read = rep(c("first", "second"), 4))
}





print("Constructing transcript ranges of all reads")

gtf <- import(GTF)

iso_results <- read.table(SIM_ISOFORMS_RESULTS, header=TRUE)
iso_results <- cbind(iso_results,sid=1:nrow(iso_results))
## sid = represent which transcript this read is simulated from.

## read the BAM file
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
## transcriptome ranges for read 1 and read 2
tr_ranges1 <- IRanges(start=read_tr_range$start1, end=read_tr_range$end1,
                      names=read_tr_range$transcript_id)
mcols(tr_ranges1)$read_id <- read_tr_range$read_id
tr_ranges2 <- IRanges(start=read_tr_range$start2, end=read_tr_range$end2,
                      names=read_tr_range$transcript_id)
mcols(tr_ranges2)$read_id <- read_tr_range$read_id

# ### generate and load ensemblDB from the GTF file
# if(!file.exists(ENSEMBL_DB_SQLITE)){
#   edb <- ensDbFromGtf(gtf = GTF, organism="Homo_sapiens", genomeVersion="GRCh37",
#                       version="85", outfile = ENSEMBL_DB_SQLITE)
# } else{
#   edb <- EnsDb(ENSEMBL_DB_SQLITE)
# }


print("Mapping transcriptomic to genomic coordinates")
## Map the transcriptomic coordinates to genomic coordinates -------------------

## Alternative to "transcriptToGenome":

## mapFromTranscripts from the GenomicFeatures packages:
## we split the reads in seperate ranges based on their splice junctions

## we identify all reads that overlap exon-exon boundaries on the transcript
## to do that, we identify the exon-exon boundaries in transcriptomic coordinates
## we take the exon-exon boundary position and map it to genomic coordinates
## we extract the SJ from the "N" positions in the CIGAR string
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

## keep all reads that overlap an exon-exon junction
# r1_truth_sj <- subsetByOverlaps(split(tr_ranges1, names(tr_ranges1)), tr_eej)

## get the eej per read, together with the read_id
tr_ranges1_list <- split(tr_ranges1, names(tr_ranges1))
tr_ranges2_list <- split(tr_ranges2, names(tr_ranges2))
r1_hits <- findOverlaps(tr_ranges1_list, tr_eej)
r2_hits <- findOverlaps(tr_ranges2_list, tr_eej)
## remove all transcripts without a read that overlaps the eej
r1_hits <- r1_hits[lengths(r1_hits)>0]
r2_hits <- r2_hits[lengths(r2_hits)>0]
## GRanges of all exon exon junctions that overlap with each read
r1_truth_eej <- lapply(seq_along(r1_hits), function(i) 
  get_eej_df(names(r1_hits)[i], r1_hits[[i]], tr_ranges1_list) )
r1_truth_eej <- GRanges(rbindlist(r1_truth_eej))

r2_truth_eej <- lapply(seq_along(r2_hits), function(i) 
  get_eej_df(names(r2_hits)[i], r2_hits[[i]], tr_ranges2_list) )
r2_truth_eej <- GRanges(rbindlist(r2_truth_eej))


## genomic coordinates of the exon exon junctions (including last and first base of exon)
##  XXXXX---XXXXXX  annotation
##      xxxxx       genomic eej coordinates
read_id1 <- mcols(r1_truth_eej)$read_id
r1_truth_eej <- mapFromTranscripts(r1_truth_eej, tr_granges)
read_id2 <- mcols(r2_truth_eej)$read_id
r2_truth_eej <- mapFromTranscripts(r2_truth_eej, tr_granges)
## add the read name
mcols(r1_truth_eej)$read_id <- read_id1[mcols(r1_truth_eej)$xHits]
mcols(r2_truth_eej)$read_id <- read_id2[mcols(r2_truth_eej)$xHits]

## remove the last and firts base of the touching exons from the range, we only 
## want the sj
r1_truth_eej <- narrow(r1_truth_eej, start=2, end=-2)
r2_truth_eej <- narrow(r2_truth_eej, start=2, end=-2)



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

## convert to GRanges
s1_range <- unlist(s1_range)
s2_range <- unlist(s2_range)

s1_idx <- match(names(s1_range), names(s1))
s1_range <- GRanges(seqnames(s1[s1_idx]), 
                     ranges=s1_range, 
                     strand(strand(s1[s1_idx])))

s2_idx <- match(names(s2_range), names(s2))
s2_range <- GRanges(seqnames(s2[s2_idx]), 
                    ranges=s2_range, 
                    strand(strand(s2[s2_idx])))

## The strand of the pair is the strand of its last alignment --> revert first
strand(s1_range) <- strand(GenomicAlignments::second(aln)[match(names(s1_range), 
                                                                names(aln))])

##----------------------------------------------------------------------------
## compare the read gaps to the true eej per read
print("Evaluating all reads")

## split the truth and the gaps per read 
## check if all ranges in query overlap with a range in truth
s1_range_list <- split(s1_range, names(s1_range))
r1_truth_eej_list <- split(r1_truth_eej, r1_truth_eej$read_id)
s2_range_list <- split(s2_range, names(s2_range))
r2_truth_eej_list <- split(r2_truth_eej, r2_truth_eej$read_id)


# res <- evaluat_read_sj_slow(s1_range_list, s2_range_list, 
                            # r1_truth_eej_list, r2_truth_eej_list, 
                            # aln)
res <- evaluat_read_sj(s1_range, s2_range, 
                       r1_truth_eej, r2_truth_eej, 
                       aln)
write.table(res, file = file.path(paste0(OUTPREFIX, "_evaluation_SJ_all.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE)




##----------------------------------------------------------------------------
## We want all reads that were simulated from one of the exons that were
## removed from the gtf annotation.
r_gtf <- import(REMOVED_GTF)

## only keep the read pairs where any of the single reads overlaps with the 
## location of the removed exons
print("Only reads from the removed exons")

read_id1 <- mcols(tr_ranges1)$read_id
read_id2 <- mcols(tr_ranges2)$read_id

r1_truth <- mapFromTranscripts(GRanges(names(tr_ranges1), tr_ranges1), tr_granges)
r2_truth <- mapFromTranscripts(GRanges(names(tr_ranges2), tr_ranges2), tr_granges)
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

## convert to GRanges
s1_range <- unlist(s1_range)
s2_range <- unlist(s2_range)

s1_idx <- match(names(s1_range), names(s1))
s1_range <- GRanges(seqnames(s1[s1_idx]), 
                    ranges=s1_range, 
                    strand(strand(s1[s1_idx])))

s2_idx <- match(names(s2_range), names(s2))
s2_range <- GRanges(seqnames(s2[s2_idx]), 
                    ranges=s2_range, 
                    strand(strand(s2[s2_idx])))

## The strand of the pair is the strand of its last alignment --> revert first
strand(s1_range) <- strand(GenomicAlignments::second(aln)[match(names(s1_range), 
                                                                names(aln))])

##----------------------------------------------------------------------------
## compare the read gaps to the true eej per read
print("Removed: Evaluating all reads")
## split the truth and the gaps per read 
## check if all ranges in query overlap with a range in truth
s1_range_list <- split(s1_range, names(s1_range))
r1_truth_eej <- r1_truth_eej[mcols(r1_truth_eej)$read_id %in% r_removed]
r1_truth_eej_list <- split(r1_truth_eej, r1_truth_eej$read_id)

s2_range_list <- split(s2_range, names(s2_range))
r2_truth_eej <- r2_truth_eej[mcols(r2_truth_eej)$read_id %in% r_removed]
r2_truth_eej_list <- split(r2_truth_eej, r2_truth_eej$read_id)

# res <- evaluat_read_sj_slow(s1_range_list, s2_range_list, 
#                        r1_truth_eej_list, r2_truth_eej_list, 
#                        aln_removed)

res <- evaluat_read_sj(s1_range, s2_range, 
                       r1_truth_eej, r2_truth_eej, 
                       aln_removed)
write.table(res, file = file.path(paste0(OUTPREFIX, "_evaluation_SJ_overl_removed_exons.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE)

##----------------------------------------------------------------------------











## TODO: speed up
## or use setdiff??
# setdiff(s1_range_list[-fp1], r1_truth_eej_list[names(s1_range_list[-fp1])])
## or mapply
# mapply(function(x,y) x %in% y, q1, t1) 
## or simply compare the gaps in each read to the list of all eej?
# s1_range_list[-fp1] %in% r1_truth_eej
## compare all gaps to all eej and count the read names with FALSE gaps


## we only count the read gaps that cannot be matched with a true sj
## wrong gaps:
# f1 <- unique(names(s1_range)[which(!names(s1_range) %in% mcols(r1_truth_eej)$read_id)])
# f2 <- unique(names(s2_range)[which(!names(s2_range) %in% mcols(r2_truth_eej)$read_id)])
# 


# s1_true_gap <- overlapsAny(s1_range, r1_truth_eej, type="equal")
# s2_true_gap <- overlapsAny(s2_range, r2_truth_eej, type="equal")
# fp1 <- length(unique(names(s1_range)[ names(s1_range) %in% 
#                                       unique(c(f1, names(s1_range)[!s1_true_gap]))]))
# fp2 <- length(unique(names(s2_range)[ names(s2_range) %in% 
#                                       unique(c(f2, names(s2_range)[!s2_true_gap]))]))



# ## TODO: speed this up?
# comp1 <- lapply(seq_along(r1_mapped), function(i) 
#   overlapsAny(r1_mapped[[i]], r1_truth[[i]], type="equal"))
# 
# ## get the list of all wrong reads
# comp2 <- setdiff(r1_mapped, r1_truth)
# 
# ## Do both methods give the same list of reads with wrong mappings?
# comp1_wrong <- names(r1_mapped)[ !sapply(comp1, all) ]
# comp2_wrong <- unique(names(unlist(comp2)))
# setequal(comp1_wrong, comp2_wrong)



## conda env for R-3.5.1 is  /home/Shared/kathi/microexon_pipeline/.snakemake/conda/a0a0dd0e
# R_LIBS=/home/Shared/Rlib/release-3.5-lib/ /usr/local/R/R-3.4.0/bin/R
