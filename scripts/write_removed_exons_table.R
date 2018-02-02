### This script constructs a summary data frame that contains information about all the novel exons in the simulated data.
## This can be used to evaluate the performance of the microexons prediction tools


library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)


OUTFILE <- snakemake@output
REMOVED_GTF <- snakemake@input[["removed_gtf"]] 
TF <- snakemake@input[["truth"]]
REDUCED_GTF <- snakemake@input[["reduced_gtf"]] 



OUTFILE <- "simulation/analysis/removed_exon_truth/removed_me_summary_table.txt"
REMOVED_GTF <- "simulation/reduced_GTF/removed_me.gtf"
TF <- "simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt"
REDUCED_GTF <- "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf" 





# rule removed_exons_table:
#     input: simulation/reduced_GTF/removed_me.gtf, simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt, simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf
#     output: simulation/analysis/removed_exon_truth/removed_me_summary_table.txt
#     jobid: 7
#     reason: Missing output files: simulation/analysis/removed_exon_truth/removed_me_summary_table.txt
#     wildcards: removed_exon=me


# OUTDIR <- "test_reduce_gtf"
# REMOVED_GTF <- file.path(OUTDIR, "removed_exons.gtf")
# GTF <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/include_mic_full_coordinate_gtf/Homo_sapiens.GRCh37.85_chr19_22_me_included_all_trans_quoted.gtf"
# tf <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/exons_truth/exons_truth_all_exons.txt"
# tf <- "simulation/analysis/GRCh37.85_chr19_22_exon_truth.txt"
# reduced_gtf <- import("test_reduce_gtf/GRCh37.85_chr19_22_reduced_exon.gtf")


reduced_gtf <- import(REDUCED_GTF)
novel_exons <- import(REMOVED_GTF)
truth <- read.table(TF, head= TRUE)



me_count <- length(novel_exons)
df <- data.frame(seqnames=seqnames(novel_exons), sj_up= rep(-1, me_count), start=start(novel_exons), end=end(novel_exons), sj_down=rep(-1, me_count), length=width(novel_exons), count=rep(NA, me_count), unique_exon = rep(NA, me_count), shared_exon_start = rep(0, me_count), shared_exon_end = rep(0, me_count), exon_overlap = rep(0, me_count))


## all exons including the removed ones
mcols(novel_exons)$count_reads <- NULL
reduced_exons <- reduced_gtf[mcols(reduced_gtf)$type == "exon"]
mcols(reduced_exons)[,c("protein_id", "protein_version")] <- NULL
exons_gtf <- c(reduced_exons, novel_exons)



# exons_gtf <- reduced_gtf[mcols(reduced_gtf)$type == "exon"]


### get the end and start coordinates of the upstream and downstream exons
## get the transcript annotation
## take the exons -1 and +1 exon_nr
## get their end and start coordinates
## if the exon is the first or last in the transcript
## TODO: put NA in there????
get_up_downstream_exon <- function(me, exons_gtf){
	tr <- exons_gtf[mcols(exons_gtf)$transcript_id==mcols(me)$transcript_id]
	print(me)
	print(tr)
	# up <- -1
	# down <- -1
	## if the me is the first or last exon of the transcript, the up or downstream exons will be -1
	if(as.character(strand(me) )== "+"){
		# print(end(tr[ mcols(tr)$exon_number==as.numeric(mcols(me)$exon_number) - 1 ] ) == NA)
		# print(start(tr[ mcols(tr)$exon_number==as.numeric(mcols(me)$exon_number) + 1 ] ) == NA)

		up <- end(tr[ mcols(tr)$exon_number== (as.numeric(mcols(me)$exon_number) - 1) ] )
		down <- start(tr[ mcols(tr)$exon_number== (as.numeric(mcols(me)$exon_number) + 1) ] )
	} else{
		up <- end( tr[ mcols(tr)$exon_number==(as.numeric( mcols(me)$exon_number) + 1) ] )
		down <- start( tr[ mcols(tr)$exon_number==(as.numeric( mcols(me)$exon_number) - 1) ] )
	}
	print(up)
	print(down)
	print(c(up, down))
	return(c(up, down))
}



res <- lapply(novel_exons, get_up_downstream_exon, exons_gtf = exons_gtf)
print(head(res))

df$sj_up <- matrix(unlist(res), nrow=length(res), byrow=T)[,1]
df$sj_down <- matrix(unlist(res), nrow=length(res), byrow=T)[,2]

## analyse if the me is unique or if it overlaps with another exon
## if yes, what type of overlap is it?
## shared 3' or 5' ends, or just overlapping exon
get_overlap_data <- function(ind, gr, df, reduced_exons){
	me <- gr[ind,]
	olap <- findOverlaps(me, reduced_exons)
	
	if(length(olap)>0){  ## the me overlaps with another exon
		olap_exons <- reduced_exons[subjectHits(olap)]
		# nr_olaps <- length(olap_exons)
		# nr_within <- length(findOverlaps(me, olap_exons, type = "within"))
		nr_shared_start <- length(findOverlaps(me, olap_exons, type = "start"))
		nr_shared_end <- length(findOverlaps(me, olap_exons, type = "end"))
		# nr_equal <- length(findOverlaps(me, olap_exons, type = "equal"))
		# print(paste0("within: ", as.character(nr_within), ", equal: ", as.character(nr_equal), ", shared start: ", as.character(nr_shared_start), ", shared end: ", as.character(nr_shared_end), ", total: ", nr_olaps) )
		df[ind, c("unique_exon", "shared_exon_start", "shared_exon_end", "exon_overlap")] <- c(FALSE, nr_shared_start, nr_shared_end, length(olap_exons)-nr_shared_start-nr_shared_end)

	}else{
		df[ind, "unique_exon"] <- TRUE
	}
	df 
}

reduced_exons <- unique(reduced_exons)
for(i in seq_along(novel_exons)){
	df <- get_overlap_data(i, gr = novel_exons, df = df, reduced_exons = reduced_exons)
}
df$unique_exon <- as.logical(df$unique_exon)

### get the true counts for all exons
tkey <- with(truth, paste0(seqnames,":",start,"-",end))
dfkey <- paste0(df$seqnames, ":", df$start, "-", df$end)
m <- match(dfkey, tkey)
df$count <- truth[m, "count_reads"]
df <- unique(df)

write.table(df, file = OUTFILE, sep = "\t", row.names = FALSE, quote = FALSE)




#### TODO !!!!! run on BAM files#####
### get the mapped counts for the upstream and downstream junctions --> can we actually predict this exon?
# BAM <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/data/STAR_2pass_outSJfilterOverhangMin6/chr19_22_SRR3192428_pass2_s.bam"


## TODO: only on the unique novel exons?? otherwise they do not match the exons in df anymore


# param <- ScanBamParam(which = novel_exons, what = c("qname"))
# reads <- readGAlignmentPairs(BAM, index = BAM, with.which_label=TRUE, param=param) ## all reads that overlap with the me regions
# ## add the me identifier to each read pair

# mcols(reads)$which_label <- mcols(first(reads))$which_label

# junc <- junctions(reads, use.mcols=TRUE) ## all junctions from the reads

# ## split the junctions according to the me
# junc <- split(junc, mcols(junc)$which_label)
# junc <- lapply(junc, unlist)
# ## the junctions start and end positions are the first and last nt in the intron (excluding the start and end!)

# ## This functions counts the number of junctions that splice to the start of the me
# count_me_start_junctions <- function(me_start, junction){
# 	sum(end(junction)+1 == me_start)
# }

# ## This functions counts the number of junctions that splice from the end of the me
# count_me_end_junctions <- function(me_end, junction){
# 	sum(start(junction)-1 == me_end)
# }


# ## for each of the me, count the junctions spanning to and from it
# df$junc_up_count <- sapply(seq_along(junc), function(x) count_me_start_junctions(df$start[x], junc[[x]]))
# df$junc_end_count <- sapply(seq_along(junc), function(x) count_me_end_junctions(df$end[x], junc[[x]]))



# ### --> in most cases, we have junctions if the exon is expressed
# ## but there are a few cases, where we have junctions, even though the true counts are 0???? --> check if this is because of overlapping exons
# ## or we don't have any junctions, even though the me is strongly expressed --> are these the very small exons???


# write.table(df, file = "../script_output/novel_me_df/novel_me_df.txt", sep = "\t", row.names = FALSE)


