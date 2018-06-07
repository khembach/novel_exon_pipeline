## This script takes the gffcompare output (from the comparison of a predicted GTF file and a reference GTF file) as input and extracts a list of all transcript IDs that are
# - correctly predicted
# - wrongly predicted
# - not predicted (missing)

library(rtracklayer)
library(data.table)
library(dplyr)
library(ggplot2)

## For the gffcompare output files and class codes see: https://ccb.jhu.edu/software/stringtie/gffcompare.shtml

# stringtie_GTF <- "/home/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/predictions/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6_stringtie.gtf"
# my_GTF <- "/home/Shared/kathi/microexon_pipeline/simulation/reduced_GTF_with_predicted_exons/me_exon/GRCh37.85_chr19_22_novel_exons_outSJfilterOverhangMin6.gtf"
# GTF <- "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"

GTF <- "/Volumes/Shared/kathi/microexon_pipeline/annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"

## the gffcompare gtf output with all transcripts
# stringtie_gffcompare_gtf <- "simulation/analysis/stringtie/gffcompare/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6/stringtie.annotated.gtf"
# prediction_gffcompared_gtf <- "simulation/analysis/gffcompare/me_exon/outSJfilterOverhangMin6/prediction.annotated.gtf"
stringtie_gffcompare_gtf <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/gffcompare/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6/stringtie.annotated.gtf"
prediction_gffcompared_gtf <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/gffcompare/me_exon/outSJfilterOverhangMin6/prediction.annotated.gtf"


## the Salmon transcript quantifications
# stringtie_qf <- "simulation/analysis/stringtie/Salmon/quantification/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6/quant.sf"
# prediction_qf <- "/home/Shared/kathi/microexon_pipeline/simulation/quantification/Salmon/me_exon/outSJfilterOverhangMin6/quant.sf"
# tpm_truth_file <- "simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results"
stringtie_qf <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/Salmon/quantification/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6/quant.sf"
prediction_qf <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/quantification/Salmon/me_exon/outSJfilterOverhangMin6/quant.sf"
tpm_truth_file <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/simulated_data/simulated_reads_chr19_22.sim.isoforms.results"


OUTDIR <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/salmon_transcript_TPM_vs_truth"


gtf <- import(GTF)
s_gtf <- import(stringtie_gffcompare_gtf)
p_gtf <- import(prediction_gffcompared_gtf)

s_gtf <- s_gtf[mcols(s_gtf)$type == "transcript"] ## we only need the transcripts
p_gtf <- p_gtf[mcols(p_gtf)$type == "transcript"] ## we only need the transcripts


### the class_code column tells us the relationship between the transcripts and the reference transcripts
## "=" indicates a complete match of the intron chain (--> the transcript are identical apart from the start of the first exon respectively the end of the last exon

## get a data frame with the transcript ids of the prediction and the corresponding matched reference id
get_matched_ids <- function(gtf){
  matched_tr <- gtf[mcols(gtf)$class_code == "="]
  return( data.frame(transcript_id = mcols(matched_tr)$transcript_id,
                     reference_id = mcols(matched_tr)$cmp_ref) )
}

## get the transcript ids of all wrong predictions
get_wrong_ids <- function(gtf){
  wrong_tr <- gtf[mcols(gtf)$class_code != "="]
  return( data.frame(transcript_id= mcols(wrong_tr)$transcript_id,
                     reference_id =mcols(wrong_tr)$cmp_ref) )
}

## get all reference ids that are missing in the prediction
get_missing_ids <- function(gtf, matched_tr_ids){
  gtf_tr <- gtf[mcols(gtf)$type == "transcript"]
  return( mcols(gtf_tr[ ! mcols(gtf_tr)$transcript_id %in%
                          matched_tr_ids$reference_id
                        ])$transcript_id )
}


## all correctly assembled transcripts:
s_matched_tr_ids <- get_matched_ids(s_gtf)
p_matched_tr_ids <- get_matched_ids(p_gtf)

## all wrongly assembled transcripts:
s_wrong_tr_ids <- get_wrong_ids(s_gtf)
p_wrong_tr_ids <- get_wrong_ids(p_gtf)

## all missing transcripts (they are in the reference, but not predicted)
s_missing_tr_ids <- get_missing_ids(gtf, s_matched_tr_ids)
p_missing_tr_ids <- get_missing_ids(gtf, p_matched_tr_ids)

#### Read the Salmon quantifications
s_quant <- fread(stringtie_qf)
p_quant <- fread(prediction_qf)
tpm_truth <- fread(tpm_truth_file)

## reference annotation transcript_id gene_id mapping
gtf_tr <- gtf[mcols(gtf)$type == "transcript"]
ref_tr_gene_id <- data.frame( mcols(gtf_tr)[,c("transcript_id", "gene_id")] )

## match the stringtie transcripts to the reference genes via the reference transcripts
s_tr_ids <- data.frame(transcript_id= mcols(s_gtf)$transcript_id,
                       reference_id =mcols(s_gtf)$cmp_ref)
## add the reference transcript id
s_quant <- s_quant %>% merge(s_tr_ids, by.x = "Name", by.y = "transcript_id")
## add the reference gene id
s_quant <- s_quant %>% merge(ref_tr_gene_id, by.x = "reference_id", by.y = "transcript_id")
p_quant <- p_quant %>% merge(ref_tr_gene_id, by.x = "Name", by.y = "transcript_id")

## Summarise the TPM counts for each gene
gene_truth <- tpm_truth %>% group_by(gene_id) %>% summarise(TPM_sum = sum(TPM))
s_gene_TPM <- s_quant %>% group_by(gene_id) %>% summarise(TPM_sum = sum(TPM))
p_gene_TPM <- p_quant %>% group_by(gene_id) %>% summarise(TPM_sum = sum(TPM))

## summarise the TPM for the TP, FP and FN transcripts
s_TP_gene_TPM <- s_quant %>%
  filter(Name %in% s_matched_tr_ids$transcript_id) %>%
  group_by(gene_id) %>%
  summarise(TPM_sum = sum(TPM))

s_FP_gene_TPM <- s_quant %>%
  filter(Name %in% s_wrong_tr_ids$transcript_id) %>%
  group_by(gene_id) %>%
  summarise(TPM_sum = sum(TPM))

p_FP_gene_TPM <- p_quant %>%
  filter(Name %in% p_wrong_tr_ids$transcript_id) %>%
  group_by(gene_id) %>%
  summarise(TPM_sum = sum(TPM))

## true tpm of the transcripts that were not predicted
s_FN_true_gene_TPM <- tpm_truth %>%
  filter(transcript_id %in% s_missing_tr_ids) %>%
  group_by(gene_id) %>%
  summarise(TPM_sum = sum(TPM))

s_TP_true_gene_TPM <- tpm_truth %>%
  filter(transcript_id %in% s_matched_tr_ids$reference_id) %>%
  group_by(gene_id) %>%
  summarise(TPM_sum = sum(TPM))

s_FP_true_gene_TPM <- tpm_truth %>%
  filter(transcript_id %in% s_wrong_tr_ids$reference_id) %>%
  group_by(gene_id) %>%
  summarise(TPM_sum = sum(TPM))

p_FN_true_gene_TPM <- tpm_truth %>%
  filter(transcript_id %in% p_missing_tr_ids) %>%
  group_by(gene_id) %>%
  summarise(TPM_sum = sum(TPM))


## transcript TPM sum of the TP (some tr are matched with the same reference)
s_TP_tr_TPM <-  s_quant %>%
  filter(Name %in% s_matched_tr_ids$transcript_id) %>%
  group_by(reference_id) %>%
  summarise(TPM_sum = sum(TPM)) %>%
  dplyr::rename(transcript_id = reference_id)

p_TP_tr_TPM <-  p_quant %>%
  merge(p_matched_tr_ids, by.x = "Name", by.y = "transcript_id") %>%
  group_by(reference_id)  %>%
  summarise(TPM_sum = sum(TPM)) %>%
  dplyr::rename(transcript_id = reference_id)


## =========== Plotting ================= ##
### Plot the true gene TPM and the stringtie TPM
gene_s_truth <- gene_truth %>% merge(s_gene_TPM, by = "gene_id", suffixes = c("_truth", "_stringtie"), all = TRUE)
gene_s_truth$TPM_sum_stringtie[ is.na(gene_s_truth$TPM_sum_stringtie) ] <- 0

p <- ggplot(gene_s_truth, aes(x = TPM_sum_stringtie, y = TPM_sum_truth)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(trans="log1p", limits=c(0,max(gene_s_truth$TPM_sum_stringtie))) +
  scale_y_continuous(trans="log1p", limits=c(0,max(gene_s_truth$TPM_sum_truth))) +
  ggtitle("Sum of transcript TPMs per gene")
ggsave(p, filename = paste0(OUTDIR, "/stringtie_vs_truth_gene_TPM.png"))

## true gene TPM and my predicted TPM
gene_p_truth <- gene_truth %>% merge(p_gene_TPM, by = "gene_id", suffixes = c("_truth", "_prediction"), all = TRUE)
gene_p_truth$TPM_sum_prediction[ is.na(gene_p_truth$TPM_sum_prediction) ] <- 0

p <- ggplot(gene_p_truth, aes(x = TPM_sum_prediction, y = TPM_sum_truth)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(trans="log1p", limits=c(0,max(gene_p_truth$TPM_sum_prediction))) +
  scale_y_continuous(trans="log1p", limits=c(0,max(gene_p_truth$TPM_sum_truth))) +
  ggtitle("Sum of transcript TPMs per gene")
ggsave(p, filename = paste0(OUTDIR, "/prediction_vs_truth_gene_TPM.png"))


## true transcript TPM and stringtie TPM of the TP transcripts
tr_s_truth <- s_TP_tr_TPM %>%
  merge( dplyr::select(tpm_truth, transcript_id, TPM), all.x=TRUE) %>%
  dplyr::rename(TPM_stringtie = TPM_sum, TPM_truth = TPM)

p <- ggplot(tr_s_truth, aes(x = TPM_stringtie, y = TPM_truth)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(trans="log1p", limits=c(0,max(tr_s_truth$TPM_stringtie))) +
  scale_y_continuous(trans="log1p", limits=c(0,max(tr_s_truth$TPM_truth))) +
  ggtitle("TPMs of TP transcripts")
ggsave(p, filename = paste0(OUTDIR, "/TP_stringtie_vs_truth_transcript_TPM.png"))

## true transcript TPM and my predicted TPM of the TP transcripts
tr_p_truth <- p_TP_tr_TPM %>%
  merge( dplyr::select(tpm_truth, transcript_id, TPM), all.x=TRUE) %>%
  dplyr::rename(TPM_prediction = TPM_sum, TPM_truth = TPM)

p <- ggplot(tr_p_truth, aes(x = TPM_prediction, y = TPM_truth)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(trans="log1p", limits=c(0,max(tr_p_truth$TPM_prediction))) +
  scale_y_continuous(trans="log1p", limits=c(0,max(tr_p_truth$TPM_truth))) +
  ggtitle("TPMs of TP transcripts")
ggsave(p, filename = paste0(OUTDIR, "/TP_prediction_vs_truth_transcript_TPM.png"))


## stringtie: FN vs FP transcripts, true TPM
s_FP_FN_gene_TPM <- s_FP_gene_TPM %>%
  merge(s_FN_true_gene_TPM, by = "gene_id", suffixes = c("_FP", "_FN"), all=TRUE)
s_FP_FN_gene_TPM$TPM_sum_FP[ is.na(s_FP_FN_gene_TPM$TPM_sum_FP) ] <- 0
s_FP_FN_gene_TPM$TPM_sum_FN[ is.na(s_FP_FN_gene_TPM$TPM_sum_FN) ] <- 0

p <- ggplot(s_FP_FN_gene_TPM, aes(x = TPM_sum_FP, y = TPM_sum_FN)) +
  geom_jitter(width = 0.3, height = 0.3, alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(trans="log1p", limits=c(0,max(s_FP_FN_gene_TPM$TPM_sum_FP))) +
  scale_y_continuous(trans="log1p", limits=c(0,max(s_FP_FN_gene_TPM$TPM_sum_FN))) +
  ggtitle("Sum of transcript TPMs per gene with FP or FN transcripts")
ggsave(p, filename = paste0(OUTDIR, "/stringtie_FP_FN_gene_TPM.png"))


## prediction: FN vs FP transcripts, true TPM
p_FP_FN_gene_TPM <- p_FP_gene_TPM %>%
  merge(p_FN_true_gene_TPM, by = "gene_id", suffixes = c("_FP", "_FN"), all=TRUE)
p_FP_FN_gene_TPM$TPM_sum_FP[ is.na(p_FP_FN_gene_TPM$TPM_sum_FP) ] <- 0
p_FP_FN_gene_TPM$TPM_sum_FN[ is.na(p_FP_FN_gene_TPM$TPM_sum_FN) ] <- 0

p <- ggplot(p_FP_FN_gene_TPM, aes(x = TPM_sum_FP, y = TPM_sum_FN)) +
  geom_jitter(width = 0.3, height = 0.3, alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(trans="log1p", limits=c(0,max(p_FP_FN_gene_TPM$TPM_sum_FP))) +
  scale_y_continuous(trans="log1p", limits=c(0,max(p_FP_FN_gene_TPM$TPM_sum_FN))) +
  ggtitle("Sum of transcript TPMs per gene with FP or FN transcripts")
ggsave(p, filename = paste0(OUTDIR, "/prediction_FP_FN_gene_TPM.png"))


##  same scale for stringtie and predictions
max_x <- max(gene_s_truth$TPM_sum_FP, gene_s_truth$TPM_sum_FP)

\
## scatter plots similar to the plotSmear from edgeR: only points that are 0 in one of the two axis are plotted as a "smear" --> with a jitter

smearWidth <- 0.1
x <- p_FP_FN_gene_TPM$TPM_sum_FP
y <- p_FP_FN_gene_TPM$TPM_sum_FN

# x <- log(x)
# y <- log(y)
w <- v <- rep(FALSE, length(x))
# w <- is.infinite(x)
# v <- is.infinite(y)
w <- x == 0
v <- y == 0
if (any(w)){
  x[w] <- min(x[!w]) - runif(sum(w), min = 0, max = smearWidth)
  x <- x - min(x) + 0.01
}
if (any(v)){
  y[v] <- min(y[!v]) - runif(sum(v), min = 0, max = smearWidth)
  y <- y - min(y) + 0.01
}

allCol <- "black"
lowCol <- "orange"
col <- rep(allCol, length(x))
if (any(w) | any(v))
  col[w | v] <- lowCol

plot(x, y, col=col, log="xy")




### same plot with ggplot2:
s_FP_FN_gene_TPM <- s_FP_gene_TPM %>%
  merge(p_FN_true_gene_TPM, by = "gene_id", suffixes = c("_FP", "_FN"), all=TRUE)
s_FP_FN_gene_TPM$TPM_sum_FP[ is.na(s_FP_FN_gene_TPM$TPM_sum_FP) ] <- 0
s_FP_FN_gene_TPM$TPM_sum_FN[ is.na(s_FP_FN_gene_TPM$TPM_sum_FN) ] <- 0

p_FP_FN_gene_TPM <- p_FP_gene_TPM %>%
  merge(p_FN_true_gene_TPM, by = "gene_id", suffixes = c("_FP", "_FN"), all=TRUE)
p_FP_FN_gene_TPM$TPM_sum_FP[ is.na(p_FP_FN_gene_TPM$TPM_sum_FP) ] <- 0
p_FP_FN_gene_TPM$TPM_sum_FN[ is.na(p_FP_FN_gene_TPM$TPM_sum_FN) ] <- 0


max_x <- max(s_FP_FN_gene_TPM$TPM_sum_FP, p_FP_FN_gene_TPM$TPM_sum_FP) + 1
max_y <- max(s_FP_FN_gene_TPM$TPM_sum_FN, p_FP_FN_gene_TPM$TPM_sum_FN) + 1

## stringtie:
smearWidth <- 0.4
w <- v <- rep(FALSE, length(s_FP_FN_gene_TPM$TPM_sum_FP))
w <- s_FP_FN_gene_TPM$TPM_sum_FP == 0
v <- s_FP_FN_gene_TPM$TPM_sum_FN == 0

if (any(w)){
  s_FP_FN_gene_TPM$TPM_sum_FP[w] <- min(s_FP_FN_gene_TPM$TPM_sum_FP[!w]) - runif(sum(w), min = 0, max = smearWidth)
  s_FP_FN_gene_TPM$TPM_sum_FP <- s_FP_FN_gene_TPM$TPM_sum_FP - min(s_FP_FN_gene_TPM$TPM_sum_FP) ## shift all values back so the minimum is 0
}
if (any(v)){
  s_FP_FN_gene_TPM$TPM_sum_FN[v] <- min(s_FP_FN_gene_TPM$TPM_sum_FN[!v]) - runif(sum(v), min = 0, max = smearWidth)
  s_FP_FN_gene_TPM$TPM_sum_FN <- s_FP_FN_gene_TPM$TPM_sum_FN - min(s_FP_FN_gene_TPM$TPM_sum_FN)
}

allCol <- "TPM > 0"
lowCol <- "TPM == 0"
col <- rep(allCol, nrow(s_FP_FN_gene_TPM))
if (any(w) | any(v))
  col[w | v] <- lowCol

p <- ggplot(s_FP_FN_gene_TPM, aes(x = TPM_sum_FP, y = TPM_sum_FN, col = col)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(trans="log1p", limits=c(0,max_x)) +
  scale_y_continuous(trans="log1p", limits=c(0,max_y)) +
  ggtitle("Sum of transcript TPMs per gene with FP or FN transcripts")
ggsave(p, filename = paste0(OUTDIR, "/stringtie_FP_FN_gene_TPM_smear.png"))

## predictions:
smearWidth <- 0.4
w <- v <- rep(FALSE, length(p_FP_FN_gene_TPM$TPM_sum_FP))
w <- p_FP_FN_gene_TPM$TPM_sum_FP == 0
v <- p_FP_FN_gene_TPM$TPM_sum_FN == 0

if (any(w)){
  p_FP_FN_gene_TPM$TPM_sum_FP[w] <- min(p_FP_FN_gene_TPM$TPM_sum_FP[!w]) - runif(sum(w), min = 0, max = smearWidth)
  p_FP_FN_gene_TPM$TPM_sum_FP <- p_FP_FN_gene_TPM$TPM_sum_FP - min(p_FP_FN_gene_TPM$TPM_sum_FP) ## shift all values back so the minimum is 0
}
if (any(v)){
  p_FP_FN_gene_TPM$TPM_sum_FN[v] <- min(p_FP_FN_gene_TPM$TPM_sum_FN[!v]) - runif(sum(v), min = 0, max = smearWidth)
  p_FP_FN_gene_TPM$TPM_sum_FN <- p_FP_FN_gene_TPM$TPM_sum_FN - min(p_FP_FN_gene_TPM$TPM_sum_FN)
}

allCol <- "TPM > 0"
lowCol <- "TPM == 0"
col <- rep(allCol, nrow(p_FP_FN_gene_TPM))
if (any(w) | any(v))
  col[w | v] <- lowCol

p <- ggplot(p_FP_FN_gene_TPM, aes(x = TPM_sum_FP, y = TPM_sum_FN, col = col)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(trans="log1p", limits=c(0,max_x)) +
  scale_y_continuous(trans="log1p", limits=c(0,max_y)) +
  ggtitle("Sum of transcript TPMs per gene with FP or FN transcripts")
ggsave(p, filename = paste0(OUTDIR, "/prediction_FP_FN_gene_TPM_smear.png"))



## Why is there no point with the max x of 12356.299 in the stringtie plot????





















### example of a wrong transcript assignment:
STRG.1594.17 and STRG.1594.16 are both matched with ENST00000196548, but STRG.1594.16 should be matched with ENST00000608843 (slightly shorted transcript, same introns)

is ENST00000608843 in the gtf file?


  c("ENST00000608843", "ENST00000196548")
--> it is correct in the stringtie_GTF (only wrong after gffcompare)--> check if all the transcripts with class "=" also have a ref_id in test







## is the reference tr id the same as in the original gtf file? and just the missing NAs are filled out?
## all ids are equal except for the ones that stringtie could not match to a transcript

test <- import(stringtie_GTF)
m <- match(s_matched_tr_ids$transcript_id, mcols(test)$transcript_id)

test_ref <- mcols(test[m])$reference_id

s_matched_tr_id <- mcols(s_matched_tr)$cmp_ref

## all transcripts with no matching transcript
s_gtf[mcols(s_gtf)$transcript_id %in% mcols(test[m[which(is.na(test_ref))]])$transcript_id ]


test[mcols(test)$transcript_id %in% mcols(test[m[which(is.na(test_ref))]])$transcript_id ]



## all transcripts that were matched to a different ref tr by gffcompare
test_tr <- test[mcols(test)$type == "transcript"]
m1 <- match(mcols(s_gtf)$transcript_id, mcols(test)$transcript_id)

test_ref_id <- mcols(test[m1])$reference_id
mcols(s_gtf)$cmp_ref == test_ref_id
# --> 159 transcripts were assigned to a different reference id!
## 1004 transcripts with NA (no ref id)


## all transcripts that were assigned to a different ref id by gffcompare
s_gtf[which(mcols(s_gtf)$cmp_ref != test_ref_id )]





## solution: take the stringtie ref id, and only  use the gffcompare id for the transcripts with NA?





## or: get the list of introns per transcript --> make a string or a vector of coordinates out of it
## (intronsBytr)
## do the same for the ref annotation
## then check if there is any ref transcript with the same intron stucture?
## --> this way we do not include the start and end of terminal exons
## but we have the same problem as before, that some of the transcripts might by ambiguous...

