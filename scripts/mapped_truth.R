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
OUTFILE <- snakemake@output[["outfile"]]

# GTF <- "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
# SIM_ISOFORMS_RESULTS <- "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results"
# OUTFILE <- "simulation/mapped_truth/STAR/me_exon/default_mapped_truth.txt"
# BAM <- "simulation/mapping/STAR/me_exon/default/pass2_Aligned.out_s.bam"


mapped_offset <- function(bam, gtf, sim_isoforms_results, verbose=FALSE){
  if(verbose) print("Loading GTF and sim.isoforms.results")
  gtf <- import(GTF)

  iso_results <- read.table(SIM_ISOFORMS_RESULTS, header=TRUE)
  iso_results <- cbind(iso_results,sid=1:nrow(iso_results))
  ## sid = represent which transcript this read is simulated from.

  if(verbose) print("Reading BAM file")
  ## read the BAM file
  aln <- readGAlignmentPairs(BAM, strandMode=2, use.names=TRUE)

  if(verbose) print("Creating transcriptome ranges from read header")
  ## get the transcriptomic ranges from each simulated read
  read_tr_range <-  str_split(string = names(aln), pattern = "_", simplify = TRUE)[,c(3, 4, 5)] 
  colnames(read_tr_range) <- c("sid", "pos", "insertL")
  read_tr_range <- apply(read_tr_range,  2, as.numeric)
  read_tr_range <- as.data.frame(read_tr_range)

  read_tr_range <- read_tr_range %>% 
    dplyr::left_join(dplyr::select(iso_results, transcript_id, length, sid), by="sid")

  read_tr_range <- read_tr_range %>% 
    dplyr::mutate(start1 = length-pos-100,  ## read length 101
                  end1 = length-pos,
                  start2 = length-pos-insertL+1,
                  end2 = length-pos-insertL+101,#or +102
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

  if(verbose) print("Comparing mapped and true read locations"
  ## Compute the offset of the true and the actual mapping start/end per read
  offset_read1 <- abs(start(read1_genome_range) - start(GenomicAlignments::first(aln))) +
                  abs(end(read1_genome_range) - end(GenomicAlignments::first(aln)))
  offset_read2 <- abs(start(read2_genome_range) - start(GenomicAlignments::second(aln))) +
                  abs(end(read2_genome_range) - end(GenomicAlignments::second(aln)))
  ## offset of the read pair
  offset <- offset_read1 + offset_read2

  ## count the occurrence of a certain offset; everything with offset >=101 is binned 
  data.frame(offset=seq(0,101), count = binCounts(offset, bx=c(seq(0,101),max(offset)) ))
}


offset_tab <- mapped_offset(BAM, GTF, SIM_ISOFORMS_RESULTS, verbose = TRUE)

print("Writing results")
write.table(offset_tab, file=OUTFILE, quote=FALSE, sep="\t", row.names=FALSE)

