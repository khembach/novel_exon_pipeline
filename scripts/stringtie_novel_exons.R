## This script takes the stringtie transcript assembly as input an compares it with an annotation to find the novel exons.

GTF <- snakemake@input[["gtf"]]
STRTIE <- snakemake@input[["strtie"]]
OUTFILE <- snakemake@output[["outfile"]]

# GTF <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me_exon.gtf"
# STRTIE <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/predictions/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6_stringtie.gtf"
# OUTFILE <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/predictions/me_exon/minReadCoverage1_minIsoformAbundance0.05/novel_exons_outSJfilterOverhangMin6_stringtie.txt"



library(rtracklayer)
library(GenomicRanges)
library(dplyr)

gtf <- import(GTF)
strtie <- import(STRTIE)


# novel <- pred[which(is.na(mcols(pred)$reference_id))]
pred <- subsetByOverlaps(strtie, gtf, type = "equal", invert = TRUE) ## all novel exons that are not yet annotated
table(mcols(pred)$type) ## 410 transcripts and 571 exons

pred_exons <- pred[mcols(pred)$type =="exon"] ## all predicted novel exons
strtie_exons <- strtie[mcols(strtie)$type =="exon"]  ## all stringtie exons


## For each novel exon, get the coordinates of the preceding and subsequent exon
pred_tab <- data.frame(pred_exons, stringsAsFactors = FALSE)
pred_tab[,c("start", "end", "exon_number")] <- lapply(pred_tab[,c("start", "end", "exon_number")], as.integer)
pred_tab$cov <- as.numeric(pred_tab$cov)



exon_tab <- data.frame(strtie[mcols(strtie)$type == "exon"])
exon_tab[,c("start", "end", "exon_number")] <- lapply(exon_tab[,c("start", "end", "exon_number")], as.integer)


## get the end of the preceding exon --> lend, if the value is NA, this means that the exon is termianl and there is no neighbouring exon on this side
pred_tab <- pred_tab %>%
  mutate(exon_number_pre = exon_number - 1, exon_number_sub=exon_number + 1) %>%
  left_join( exon_tab %>% select("transcript_id", "exon_number", "end"),
             by = c("transcript_id","exon_number_pre"="exon_number") ) %>%
  rename(lend = end.y)

## get the start of the subsequent exon --> rstart
pred_tab <- pred_tab %>%
  left_join( exon_tab %>% select("transcript_id", "exon_number", "start"),
             by = c("transcript_id","exon_number_sub"="exon_number") ) %>%
  rename(rstart = start.y, start = start.x, end = end.x)


## We need a value that tells us how likely it is that this prediction is correct --> score, but this is currently not used by stringtie, all values are 1000
## We try the coverage, maybe FPKM or TPM would also work
# cov: The average per-base coverage for the transcript or exon.
# FPKM: Fragments per kilobase of transcript per million read pairs. This is the number of pairs of reads aligning to this feature, normalized by the total number of fragments sequenced (in millions) and the length of the transcript (in kilobases).
# TPM: Transcripts per million. This is the number of transcripts from this particular gene normalized first by gene length, and then by sequencing depth (in millions) in the sample. A detailed explanation and a comparison of TPM and FPKM can be found here, and TPM was defined by B. Li and C. Dewey here.
# we rename the coverage to min_reads, because then we can use the same script to plot the PR curve
results <- pred_tab %>%
  select(seqnames, lend, start, end, rstart, strand, cov) %>%
  rename(min_reads = cov)

## remove duplicate predictions: we take the maximum coverage for identical exon predictions
results <- results %>%
  group_by(seqnames, lend, start, end, rstart, strand) %>%
  summarise(min_reads = max(min_reads))

write.table(results, file=OUTFILE, row.names=FALSE, quote=FALSE, sep="\t" )




