#### This script removes some exons from a given gtf file. It makes sure that all removed exons are expressed (TPM > 0).
## It keeps track of the removed exons and the transcripts in which the exons occurred.
## The output is a reduced GTF file with missing exons (microexon, normal exons or both). Some of the removed exons are unique, and some do overlap with other exons (more complicated to discover).
## @ author: Katharina Hembach
## 25.10.2017



GTF <- snakemake@input[["gtf"]]
EXPR <- snakemake@input[["truth"]]
OUTFILE_ME <- snakemake@output[["me"]]
OUTFILE_EXON <- snakemake@output[["exon"]]
OUTFILE_ME_EXON <- snakemake@output[["me_exon"]]
# "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf
OUTDIR <- dirname(OUTFILE_ME)




GTF <- "annotation/Homo_sapiens.GRCh37.85_chr19_22.gtf"
EXPR <- "simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt"

# # OUTFILE_ME <- "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf"
# # OUTFILE_EXON <- "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_exon.gtf"
# # OUTFILE_ME_EXON <- "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me_exon.gtf"

# OUTFILE_ME <- "test_reduce_gtf/GRCh37.85_chr19_22_reduced_me.gtf"
# OUTFILE_EXON <- "test_reduce_gtf/GRCh37.85_chr19_22_reduced_exon.gtf"
# OUTFILE_ME_EXON <- "test_reduce_gtf/GRCh37.85_chr19_22_reduced_me_exon.gtf"



# OUTDIR <- "test_reduce_gtf"







library(rtracklayer)
library(ggplot2)
library(GenomicFeatures)

#################
# read the data #
#################

gtf <- import(GTF)
expr <- read.table(EXPR, header = TRUE)

# we filter out all exons that are not "normal" exons in a transcript, because we don't want short annotations such as start codons
not_exons <- gtf[mcols(gtf)$type %in% c("start_codon", "stop_codon", "Selenocysteine", "five_prime_utr", "three_prime_utr")]

exons <- subsetByOverlaps(gtf, not_exons, type = "equal", invert = TRUE)
exons <- exons[mcols(exons)$type == "exon" & mcols(exons)$gene_biotype == "protein_coding"] # only keep exons from protein coding transcripts (remove retained introns, antisense transcripts, ...)
mcols(exons)$count_reads <- expr[subjectHits(findOverlaps(exons, GRanges(expr), type = "equal")) , "count_reads" ]  ## add the simulated counts per exon


me <- exons[width(exons)<28 & width(exons)>=3,] ##845 unique microexons
me_expressed <- me[mcols(me)$count_reads >0,] 
me_expressed_unique <- unique(me_expressed) ##176 unique expressed me, count_reads>0

normal_exons <- unique(exons[width(exons)>=28,])
exons_expressed <- normal_exons[mcols(normal_exons)$count_reads >1,] ##16258 unique expressed exons (>27nt) with TPM>0
exons_expressed_unique <- unique(exons_expressed)
# only keep the exons that are in the middle of a transcript
# get all exons per transcript and remove the first and last ones
txdb <- makeTxDbFromGFF(GTF, format="gtf")
exbytx <- exonsBy(txdb, by = "tx", use.names=TRUE) 


first_exon_id <- unique(sapply(exbytx, function(x) x[1]$exon_id))
last_exon_id <- unique(sapply(exbytx, function(x) x[length(x)]$exon_id) )

firsts <- exons(txdb)[first_exon_id] # 16021
lasts <- exons(txdb)[last_exon_id] # 16047
first_last <- unique(c(firsts, lasts)) # 29864

middle_exons <- subsetByOverlaps(exons_expressed_unique, first_last, type="equal", invert=TRUE) # 8521 

## only keep the unique exons that do not overlap with another exon ?? --> not if we want to show that the Salmon derived counts are better than featureCounts


########################
# remove the me / exon #
########################
## random selection of 100 me
set.seed(42)
remove_id_me <- sample(1:length(me_expressed_unique), 100, replace=FALSE)
## remove the exons from the corresponding genes and all transcripts
to_remove_me <- me_expressed_unique[remove_id_me,]


## random selection of 100 exons
set.seed(42)
remove_id_exon <- sample(1:length(middle_exons), 100, replace=FALSE)
## remove the exons from the corresponding genes and all transcripts
to_remove_exon <- middle_exons[remove_id_exon,]


######################################
# Quality check of removed microexons #
######################################
##### Plot the exon ount distribution of all exons, all me and all me that will be removed

## get the TPM sum from all transcripts
eid <- paste0(seqnames(exons), ":", start(exons), "-", end(exons))
meid <- paste0(seqnames(me), ":", start(me), "-", end(me))
rid <- paste0(seqnames(to_remove_me), ":", start(to_remove_me), "-", end(to_remove_me))
exons_count <- lapply(split(mcols(exons)$count_reads, eid), unique)

# get the counts for each exon
df <- data.frame( id=c(names(exons_count), rid), type = c(ifelse(names(exons_count) %in% meid, "microexon", "NOTmicroexon"), rep("removed", length(rid))), read_count = c(unlist(exons_count), unlist(exons_count[match(rid, names(exons_count)) ] ))
 )

p <- ggplot(df[df$read_count >0, ], aes(x=read_count, color = type)) + geom_density() + xlim(c(0, 1000))
ggsave(file.path(OUTDIR, "density_me_read_count_seed_42.pdf"), plot = p, device = "pdf")

### Percentage of microexons that are not unique and overlapping with another exons
hits <- findOverlaps(unique(me), unique(normal_exons))
length(unique(queryHits(hits)))  ## 603 of 939 me are overlappping

hits <- findOverlaps(to_remove_me, unique(normal_exons))
length(unique(queryHits(hits)))  ## 68 of 100 me are overlappping, 32 unique me

### length distribution is also similar
summary(width(to_remove_me))
summary(width(unique(me)))


rid <- paste0(seqnames(to_remove_exon), ":", start(to_remove_exon), "-", end(to_remove_exon))

df <- data.frame( id=c(names(exons_count), rid), type = c(rep("all exons", length(exons_count)),rep("removed", length(rid))), read_count = c(unlist(exons_count), unlist(exons_count[match(rid, names(exons_count)) ] ))
 )

p <- ggplot(df[df$read_count >0, ], aes(x=read_count, color = type)) + geom_density() + xlim(c(0, 3000))
ggsave(file.path(OUTDIR, "density_normal_exons_read_count_seed_42.pdf"), plot = p, device = "pdf")

### Percentage of removed exons that are not unique and overlapping with another exons
not_remove <- subsetByOverlaps(unique(exons), to_remove_exon,  type = "equal", invert = TRUE)
hits <- findOverlaps(to_remove_exon, not_remove)
length(unique(queryHits(hits)))  ## 53 of 100 exons are overlappping
### length distribution is also similar
summary(width(to_remove_exon))
summary(width(middle_exons))



#################################
# remove the exons from the GTF #
#################################
### ME ###
hits <- findOverlaps(to_remove_me, exons, type = "equal", ignore.strand= FALSE )
to_remove_me_all <- exons[subjectHits(hits)] ## me to remove from all transcripts that contain them (157 exons, 149 transcripts)
gtf_part_me <- gtf[mcols(gtf)$transcript_id %in% mcols(to_remove_me_all)$transcript_id,] # full transcripts that contain the me that we want to remove
gtf_unchanged_me <- gtf[!mcols(gtf)$transcript_id %in% mcols(to_remove_me_all)$transcript_id,]

gtf_part_me_reduced <- subsetByOverlaps(gtf_part_me, to_remove_me_all,  type = "equal", invert = TRUE) ## keep everything that is not one of the 100 me

## ?? Do we need to do this?
## make sure the gtf part is still invalid GTF format
## check the exon_number
## make sure that all transcripts are distinct and the me was not the only difference between transcripts, if yes should we remove the transcript, or keep two annoatations for the same tr?

gtf_me_reduced <- c(gtf_unchanged_me, gtf_part_me_reduced)
# write reduced gtf to file
export(gtf_me_reduced, OUTFILE_ME)
# gtf_new <- import(OUTFILE)  ## test if it is in valid format
export(to_remove_me_all, file.path(OUTDIR, "removed_me.gtf"))


#### EXON ####
hits <- findOverlaps(to_remove_exon, exons, type = "equal", ignore.strand= FALSE )
to_remove_exon_all <- exons[subjectHits(hits)] ## me to remove from all transcripts that contain them (348 exons, 316 transcripts)
gtf_part_exon <- gtf[mcols(gtf)$transcript_id %in% mcols(to_remove_exon_all)$transcript_id,] # full transcripts that contain the me that we want to remove
gtf_unchanged_exon <- gtf[!mcols(gtf)$transcript_id %in% mcols(to_remove_exon_all)$transcript_id,]
gtf_part_exon_reduced <- subsetByOverlaps(gtf_part_exon, to_remove_exon_all,  type = "equal", invert = TRUE) #
gtf_exon_reduced <- c(gtf_unchanged_exon, gtf_part_exon_reduced)

export(gtf_exon_reduced, OUTFILE_EXON )
export(to_remove_exon_all, file.path(OUTDIR, "removed_exon.gtf"))


### ME + EXON ####
to_remove_me_exon <- c(to_remove_me, to_remove_exon)
hits <- findOverlaps(to_remove_me_exon, exons, type = "equal", ignore.strand= FALSE )
to_remove_all <- exons[subjectHits(hits)] ## me to remove from all transcripts that contain them (157 exons, 149 transcripts)
gtf_part <- gtf[mcols(gtf)$transcript_id %in% mcols(to_remove_all)$transcript_id,] # full transcripts that contain the me that we want to remove
gtf_unchanged <- gtf[!mcols(gtf)$transcript_id %in% mcols(to_remove_all)$transcript_id,]

gtf_part_reduced <- subsetByOverlaps(gtf_part, to_remove_all,  type = "equal", invert = TRUE) ## keep everything that is not one of the 100 me

gtf_me_exon_reduced <- c(gtf_unchanged, gtf_part_reduced)
export(gtf_me_exon_reduced, OUTFILE_ME_EXON)
export(to_remove_all, file.path(OUTDIR, "removed_me_exon.gtf"))



################################
# Summary of removed ME / EXON #
################################

## This function returns the two splice junctions of a exon
## NA empty string if the exon is the first of last and there is no splice junction to the left or right
## input paramters:
## index in summary data frame, exons of interest as GRange, row of summary data frame, gtf annotation, unique set of all exons in annotation
get_removed_exons_summary <- function(ind, e, df, anno, exons){
	trans <- anno[which(mcols(anno)$transcript_id == mcols(e)$transcript_id & mcols(anno)$type == "exon")]  ## all exons from the transcript
	trans <- sort(trans) # sorted with increasing coordinates (independent of transcript strand)
	me_ind <- which(start(trans) == start(e) & end(trans) == end(e))
	df$lend <- ifelse(me_ind - 1 > 0 , end(trans)[me_ind -1], NA)  ## end of left exon
	df$rstart <- ifelse(me_ind + 1 <= length(trans) , start(trans)[me_ind + 1], NA) ## start of right exons

	## number of overlapping exons (same start, end coordinates or overlappign with different ends)
	exons <- subsetByOverlaps(exons, to_remove_all,  type = "equal", invert = TRUE) # remove the exon from the annotation
	olap <- findOverlaps(e, exons)
	
	if(length(olap)>0){  ## the me overlaps with another exon
		olap_exons <- exons[subjectHits(olap)]
		nr_shared_start <- length(findOverlaps(e, olap_exons, type = "start"))
		nr_shared_end <- length(findOverlaps(e, olap_exons, type = "end"))

		df[c("unique_exon", "shared_exon_start", "shared_exon_end", "exon_overlap")] <- c(FALSE, nr_shared_start, nr_shared_end, length(olap_exons)-nr_shared_start-nr_shared_end)

	}else{
		df["unique_exon"] <- TRUE
	}	
	return(df)
}



### ME ####
## data frame with all exons that will be removed and information about the previous and sequential exon and the number of overlapping exons
me_summary <- as.data.frame(to_remove_me_all)[,c("seqnames", "start", "end", "width", "strand", "gene_id", "transcript_id", "exon_id", "count_reads")]
me_summary$lend <- NA
me_summary$rstart <- NA
me_summary$unique_exon <- NA
me_summary$shared_exon_start <- 0
me_summary$shared_exon_end <- 0
me_summary$exon_overlap <- 0

for(i in 1:nrow(me_summary)){
	me_summary[i,] <- get_removed_exons_summary(i, e = to_remove_me_all[i], df = me_summary[i,], anno = gtf, exons = unique(exons))
}
me_summary$unique_exon <- as.logical(me_summary$unique_exon)

write.table(me_summary, file = file.path(OUTDIR, "removed_microexons.txt"), sep = "\t", row.names = FALSE, quote=FALSE )
## get a list of unique splice junction combination (remove duplicate ones)
me_summary_unique <- me_summary[, !colnames(me_summary) %in% c("exon_id", "transcript_id")]
me_summary_unique <- unique(me_summary_unique)  ## 124 unique junctions
write.table(me_summary_unique, file = file.path(OUTDIR, "removed_microexons_unique.txt"), sep = "\t", row.names = FALSE, quote=FALSE )


#### EXON ####
exon_summary <- as.data.frame(to_remove_exon_all)[,c("seqnames", "start", "end", "width", "strand", "gene_id", "transcript_id", "exon_id", "count_reads")]
exon_summary$lend <- NA
exon_summary$rstart <- NA
exon_summary$unique_exon <- NA
exon_summary$shared_exon_start <- 0
exon_summary$shared_exon_end <- 0
exon_summary$exon_overlap <- 0

for(i in 1:nrow(exon_summary)){
	exon_summary[i,] <- get_removed_exons_summary(i, e = to_remove_exon_all[i], df = exon_summary[i,], anno = gtf, exons = unique(exons))
}
exon_summary$unique_exon <- as.logical(exon_summary$unique_exon)

write.table(exon_summary, file = file.path(OUTDIR, "removed_exons.txt"), sep = "\t", row.names = FALSE, quote=FALSE )
## get a list of unique splice junction combination (remove duplicate ones)
exon_summary_unique <- exon_summary[, !colnames(exon_summary) %in% c("exon_id", "transcript_id")]
exon_summary_unique <- unique(exon_summary_unique)  ## 107 unique junctions
write.table(exon_summary_unique, file = file.path(OUTDIR, "removed_exons_unique.txt"), sep = "\t", row.names = FALSE, quote=FALSE )


## ME + EXON ##
me_exon_summary <- as.data.frame(to_remove_all)[,c("seqnames", "start", "end", "width", "strand", "gene_id", "transcript_id", "exon_id", "count_reads")]
me_exon_summary$lend <- NA
me_exon_summary$rstart <- NA
me_exon_summary$unique_exon <- NA
me_exon_summary$shared_exon_start <- 0
me_exon_summary$shared_exon_end <- 0
me_exon_summary$exon_overlap <- 0

for(i in 1:nrow(me_exon_summary)){
	me_exon_summary[i,] <- get_removed_exons_summary(i, e = to_remove_all[i], df = me_exon_summary[i,], anno = gtf, exons = unique(exons))
}
me_exon_summary$unique_exon <- as.logical(me_exon_summary$unique_exon)

write.table(me_exon_summary, file = file.path(OUTDIR, "removed_microexons_exons.txt"), sep = "\t", row.names = FALSE, quote=FALSE )
## get a list of unique splice junction combination (remove duplicate ones)
me_exon_summary_unique <- me_exon_summary[, !colnames(me_exon_summary) %in% c("exon_id", "transcript_id")]
me_exon_summary_unique <- unique(me_exon_summary_unique)  ## 107 unique junctions
write.table(me_exon_summary_unique, file = file.path(OUTDIR, "removed_microexons_exons_unique.txt"), sep = "\t", row.names = FALSE, quote=FALSE )
