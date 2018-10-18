## Compare the mapped position of the simulated reads with their true location as indicated in the read header.
## Soft clipped reads are considered as wrong.

library(rtracklayer)
library(GenomicAlignments)
library(stringr)
library(dplyr)
library(GenomicFeatures)

## parameters
BAM <- snakemake@input[["bam"]]
GTF <- snakemake@input[["gtf"]]
SIM_ISOFORMS_RESULTS <- snakemake@input[["sim_iso_res"]]
OUTPREFIX <- snakemake@params[["outprefix"]]
REMOVED_GTF <- snakemake@input[["removed_gtf"]]

# GTF <- "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
# SIM_ISOFORMS_RESULTS <- "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results"
# OUTFILE <- "simulation/mapped_truth/STAR/me_exon/default_mapped_truth.txt"
# BAM <- "simulation/mapping/STAR/me_exon/default/pass2_Aligned.out_s.bam"
#  REMOVED_GTF <- "simulation/reduced_GTF/removed_me_exon.gtf"
# OUTPREFIX <- "simulation/mapped_truth/star/me_exon/default/star_"

#' Offset of aligned reads and their true location
#' 
#' Create an integer vector with the offset per read pair. The offset is 
#' the difference of read start and end of the forward and reverse reads
#' with the true genomic location of the reads.
#'
#' @param read1_truth GRanges with the true location of read1
#' @param read2_truth GRanges with the true lcoation of read2 
#' @param aln GAlignmentPairs with the aligned paired-end reads
#'
#' @return integer vector with the offset per read pair
#' @export
#'
#' @examples
compute_read_offset <- function(read1_truth, read2_truth, aln){
  offset_read1 <- abs(start(read1_truth) - start(GenomicAlignments::first(aln))) +
                  abs(end(read1_truth) - end(GenomicAlignments::first(aln)))
  offset_read2 <- abs(start(read2_truth) - start(GenomicAlignments::second(aln))) +
                  abs(end(read2_truth) - end(GenomicAlignments::second(aln)))
  ## offset of the read pair
  offset_read1 + offset_read2
}

#' Soft clipped bases per read pair
#' 
#' Compute the number of soft clipped bases per read pair.
#'
#' @param aln GAlignmentPairs with the aligned paired-end reads
#'
#' @return integer vector with the number of soft clipped bases per read pair
#' @export
#'
#' @examples
soft_clipped_per_pair <- function(aln){
  sc_read1 <- cigar( GenomicAlignments::first(aln)) %>% 
    cigarOpTable
  sc_read1 <- sc_read1[,"S"]
  
  sc_read2 <- cigar( GenomicAlignments::second(aln)) %>% 
    cigarOpTable
  sc_read2 <- sc_read2[,"S"]
  
  sc_read1 + sc_read2
}


#' Mapped truth of BAM file
#'
#' Computes the sum of offsets between the genomic locations of a read-pair in 
#' the BAM file with the true simulated location. The offsets from all read-pairs
#' are summarized in a count table.
#' 
#' @param bam BAM file
#' @param gtf GTF file with genome annotation
#' @param sim_isoforms_results SAMPLE.sim.isoforms.results output from RSEM for 
#' the BAM file sample
#' @param verbose print status messages
#'
#' @return data.frame with count table of offsets
#' @export
#'
#' @examples
write_offset_tables <- function(bam, gtf_file, sim_isoforms_results, removed_gtf_file, outprefix, verbose=FALSE){
  if(verbose) print("Loading GTF and sim.isoforms.results")
  gtf <- import(gtf_file)

  iso_results <- read.table(sim_isoforms_results, header=TRUE)
  iso_results <- cbind(iso_results,sid=1:nrow(iso_results))
  ## sid = represent which transcript this read is simulated from.

  if(verbose) print("Reading BAM file")
  ## read the BAM file
  aln <- readGAlignmentPairs(bam, strandMode=2, use.names=TRUE)

  if(verbose) print("Creating transcriptome ranges from read header")
  ## get the transcriptomic ranges from each simulated read
  read_tr_range <-  str_split(string = names(aln), pattern = "_", simplify=TRUE)[,c(3, 4, 5)] 
  colnames(read_tr_range) <- c("sid", "pos", "insertL")
  read_tr_range <- apply(read_tr_range,  2, as.numeric)
  read_tr_range <- as.data.frame(read_tr_range)

  read_tr_range <- read_tr_range %>% 
    dplyr::left_join(dplyr::select(iso_results, transcript_id, length, sid), by="sid")

  read_tr_range <- read_tr_range %>% 
    dplyr::mutate(start1 = length-pos-100,  ## read length 101
                  end1 = length-pos,
                  start2 = length-pos-insertL+1,
                  end2 = length-pos-insertL+101
                  )
  ## transcriptome ranges for read 1 and read 2
  tr_ranges1 <- GRanges(read_tr_range$transcript_id, 
                        ranges=IRanges(read_tr_range$start1, read_tr_range$end1))
  tr_ranges2 <- GRanges(read_tr_range$transcript_id, 
                        ranges=IRanges(read_tr_range$start2, read_tr_range$end2))

  if(verbose) print("Making TxDb and extracting transcript annotations")
  ## get GRanges from all transcripts
  txdb <- makeTxDbFromGRanges(gtf)
  tr_granges <- exonsBy(txdb, by="tx", use.names=TRUE)

  if(verbose) print("Mapping transcript to genome ranges")
  read1_genome_range <- mapFromTranscripts(tr_ranges1, tr_granges)
  read2_genome_range <- mapFromTranscripts(tr_ranges2, tr_granges)

  ## get the number of soft clipped bases per read
  sc_pair <- soft_clipped_per_pair(aln)
  
  if(verbose) print("Filter reads from removed exons")
  ## We want all reads that were simulated from one of the exons that were
  ## removed from the gtf annotation.
  r_gtf <- import(removed_gtf_file)

  ## only keep the read pairs where any of the single reads overlaps with the location of the removed exons
  reads1_r <- subsetByOverlaps(read1_genome_range, r_gtf)
  reads2_r <- subsetByOverlaps(read2_genome_range, r_gtf)
  r_ind <- unique(c(reads1_r$xHits, reads2_r$xHits))
  aln_r <- aln[r_ind]
  read1_r_genome_range <- read1_genome_range[r_ind]
  read2_r_genome_range <- read2_genome_range[r_ind]
  
  if(verbose) print("Comparing mapped and true read locations")
  ## Compute the offset of the true and the actual mapping start/end per read
  offset <- compute_read_offset(read1_genome_range, read2_genome_range, aln)
  ## count the occurrence of a certain offset; everything with offset >=101 is binned 
  offset_tab<- data.frame(offset = seq(0,101), 
             count = binCounts(offset, bx=c(seq(0,101), max(offset))))
  write.table(offset_tab, 
              file=paste0(outprefix, "offset_counts.txt"), 
              quote=FALSE, sep="\t", row.names=FALSE)
  write.table(data.frame(offset = offset, soft_clipped = sc_pair), 
              file=paste0(outprefix, "offset_soft_clipped.txt"), 
              quote=FALSE, sep="\t", row.names=FALSE)
  ## reads from removed exons
  offset_r <- compute_read_offset(read1_r_genome_range, read2_r_genome_range, aln_r)
  offset_tab_r <- data.frame(offset_r = seq(0,101), 
                          count = binCounts(offset, bx=c(seq(0,101), max(offset))))
  write.table(offset_tab_r, 
              file=paste0(outprefix, "offset_counts_removed_exons.txt"), 
              quote=FALSE, sep="\t", row.names=FALSE)
  write.table(data.frame(offset = offset_r, soft_clipped = sc_pair[r_ind]), 
              file=paste0(outprefix, "offset_soft_clipped_removed_exons.txt"), 
              quote=FALSE, sep="\t", row.names=FALSE)
} 


write_offset_tables(BAM, GTF, SIM_ISOFORMS_RESULTS, REMOVED_GTF, OUTPREFIX, verbose=TRUE)
