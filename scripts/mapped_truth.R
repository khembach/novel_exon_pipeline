# This script computes the sids for each exon and then extracts all reads from the bam file that contain the sids. The exon read counts are them computed from the filtered bam file.

## The output file format is the following:
# nr	seqnames	start	end	length	count_reads

## only works if R is started with 'R_LIBS=/home/Shared/Rlib/release-3.2-lib/ /usr/local/R/R-3.2.2/bin/R'

library(rtracklayer)
library(GenomicAlignments)
library(Rsamtools)
library(stringr)

## parameters
BAM <- snakemake@input[["bam"]]
GTF <- snakemake@input[["gtf"]]
SIM_ISOFORMS_RESULTS <- snakemake@input[["sim_iso_res"]]
OUT_DIR <- snakemake@input[["outdir"]]



# GTF <- "/Volumes/Shared/kathi/microexon_pipeline/annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
# SIM_ISOFORMS_RESULTS <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results"
# OUT_DIR <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/mapped_truth/"



# GTF <- "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
# SIM_ISOFORMS_RESULTS <- "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results"
# OUTFILE <- "simulation/mapped_truth/STAR/me_exon/decault_mapped_truth.txt"
# BAM <- "simulation/mapping/STAR/me_exon/default/pass2_Aligned.out_s.bam"





gtf <- import(GTF)

# SIM_ISOFORMS_RESULTS 
iso.results <- read.table(SIM_ISOFORMS_RESULTS, header=TRUE)
iso.results <- cbind(iso.results,sid=1:nrow(iso.results))
## sid = represent which transcript this read is simulated from.

### GTF file with 
gtf <- gtf[mcols(gtf)$type=='exon']
mcols(gtf)$sid <- merge(as.data.frame(gtf),iso.results[,c('transcript_id','sid')],by='transcript_id',sort=F)$sid   

exons <- subset(as.data.frame(gtf),select=c(seqnames,start,end,strand, exon_id,transcript_id,sid))

# exons_range <- GRanges(seqnames="19", ranges=unique(IRanges(exons[exons$seqnames=="19",]$start, exons[exons$seqnames=="19",]$end)) )
# exons_range <- c(exons_range, GRanges(seqnames="22", ranges=unique(IRanges(exons[exons$seqnames=="22",]$start, exons[exons$seqnames=="22",]$end)) ) )
# exons_range <- sort(exons_range)
# # exons_range<- exons_range[1:2000,]   ###for testing


exons_range <- sort(unique(granges(gtf)))
exons_range$eid <- paste0(seqnames(exons_range), "_", start(exons_range), "_", end(exons_range), "_", strand(exons_range))

# exon ids and corresponding sids
eid_sid <- data.frame(eid = paste0(exons$seqnames, "_", exons$start, "_", exons$end, "_", exons$strand), sid = exons$sid)

## read BAM file
what <- c("qname")
param <- ScanBamParam(what=what)

aln <- readGAlignments(BAM, param=param)
mcols(aln)$sid <- str_split(string = mcols(aln)$qname, pattern = "_", simplify = TRUE)[,3]
# add a metadata column with the sid: the sid is part of the qname of each read


# list of sids per exon
# table with exon coordinates, eid and sid
eid_split <- split(eid_sid$sid, eid_sid$eid)
# per eid, get all reads with the sids

## Count the number of overlapping reads per exon
eid_count <- data.frame(eid = names(eid_split), read_count = NA)
for (i in 1:length(eid_split)){
	eid_count$read_count[i] <- length( subsetByOverlaps( aln[ mcols(aln)$sid %in% eid_split[[i]],],
                                                    exons_range[exons_range$eid == names(eid_split)[i]] ) )
}


results <- exons[,1:4]
results$read_count <- eid_count$read_count


write.table(results, file=OUTFILE, quote=FALSE, sep="\t", row.names=FALSE)
