## This script takes a GTF file and salmon quantificatins as input and computes the coverage of all exons.
## For each exon we compute the coverage of all transcripts in which it occurrs and sum them up to get the estimated exon counts

# MICS_INC <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/include_mic_full_coordinate_gtf/me_included.txt"

# GTF <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/include_mic_full_coordinate_gtf/Homo_sapiens.GRCh37.85_chr19_22_me_included_all_trans_quoted.gtf"
# tf <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/exons_truth/exons_truth_all_exons.txt"
# qf <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/salmon/quant.sf"

# GTF_noHeader <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/include_mic_full_coordinate_gtf/Homo_sapiens.GRCh37.85_chr19_22_me_included_all_trans_quoted.noHeader.gtf"
# SQLITE <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/include_mic_full_coordinate_gtf/Homo_sapiens.GRCh37.85_chr19_22_me_included_all_trans_quoted.gtf.sqlite"
# fldgz <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/salmon/aux_info/fld"
# SalmonOutput <- "/home/Shared/data/seq/microexon_simulation/microexon_study/new_simulation/data/"


### packages
library(rtracklayer)
library(data.table)
library(ggplot2)
library(viridis)



GTF <- snakemake@input[["gtf"]]
tf <- snakemake@input[["truth"]]
rf <- snakemake@input[["removed"]]
qf <- snakemake@input[["quant"]]
fldgz <- snakemake@input[["fldgz"]]
SalmonOutput <- snakemake@output[["outdir"]]


# GTF <- "/home/Shared/kathi/microexon_pipeline/simulation/reduced_GTF_with_predicted_exons/me_exon/GRCh37.85_chr19_22_novel_exons_outSJfilterOverhangMin6.gtf"
# tf <- "/home/Shared/kathi/microexon_pipeline/simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt"
# qf <- "/home/Shared/kathi/microexon_pipeline/simulation/quantification/Salmon/me_exon/outSJfilterOverhangMin6/quant.sf"
# fldgz <- "/home/Shared/kathi/microexon_pipeline/simulation/quantification/Salmon/me_exon/outSJfilterOverhangMin6/aux_info/fld"
# SalmonOutput <- "/home/Shared/kathi/microexon_pipeline/simulation/analysis/derived_Salmon_counts/outSJfilterOverhangMin6/"
# rf <- "/home/Shared/kathi/microexon_pipeline/simulation/reduced_GTF/removed_microexons_exons_unique_classified.txt"


# GTF <- "/home/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/predictions/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6_stringtie.gtf"
# tf <- "/home/Shared/kathi/microexon_pipeline/simulation/analysis/GRCh37.85_chr19_22_all_exon_truth.txt"
# qf <- "/home/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/Salmon/quantification/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6/quant.sf"
# fldgz <- "/home/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/Salmon/quantification/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6/aux_info/fld.gz"
# SalmonOutput <- "/home/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/derived_Salmon_counts/me_exon/minReadCoverage1_minIsoformAbundance0.05/outSJfilterOverhangMin6/test/"
# rf <- "/home/Shared/kathi/microexon_pipeline/simulation/reduced_GTF/removed_me_exon_unique_classified.txt"



## parse the gtf file to get a list of all exons and the corresponding transcripts
gtf <- import(GTF)    
gtf <- gtf[mcols(gtf)$type=="exon"]                          # gtf with only exon types
ex <- paste0(seqnames(gtf),"_",start(gtf),":",end(gtf))      # exon locations pasted
# trans <- split(mcols(gtf)$transcript_id, ex)                 # transcripts grouped by the exon location

## read the truth and match the exon from the transcript list to the exons in the truth table 
truth <- fread(tf)
## read the Salmon quantifications
quant <- fread(qf)
# u <- unique(unlist(trans))
# quant <- quant[quant$Name %in% u,]  ## all transcripts that contain the exons


## get the locations of the novel microexons and exons
removed <- fread(rf)
me_id <- removed$type =="me"
removed_me_loc <- unique( paste0(removed$seqnames[me_id], "_", removed$start[me_id], ":", removed$end[me_id]) )
exon_id <- removed$type =="exon"
removed_exon_loc <- unique( paste0(removed$seqnames[exon_id], "_", removed$start[exon_id], ":", removed$end[exon_id]) )


# novel_exons <- gtf[grepl("_exon_", mcols(gtf)$exon_id)]
# novel_me <- gtf[grepl("_me_", mcols(gtf)$exon_id)]


### Compute the mean fragment length from the Salmon auxilliary files
## fld.gz: This file contains an approximation of the observed fragment length distribution. It is a gzipped, binary file containing integer counts. The number of (signed, 32-bit) integers (with machine-native endianness) is equal to the number of bins in the fragment length distribution (1,001 by default — for fragments ranging in length from 0 to 1,000 nucleotides).
system(paste0("gunzip -c ",fldgz, " > ", file.path(SalmonOutput, "fld" )) )
fld <- file.path(SalmonOutput, "fld" )

flBins <- readBin(fld, integer(), n=1001, signed=TRUE)
bins <- 0:1000
fragmentLength <-bins[which.max(flBins)]
print(paste0("Fragment length: ", fragmentLength))


### compute the start and end position of all exons in the transcripts and add them to the gtf file
## we assume that all exons are ordered ascending according to their location on the genome (independent of strand

## we assume that all exons in a "+"-strand transcript are ordered ascending, and decreasing on the "-" strand
## --> we split the gtf and order the two strands seperately, then we join it together again
gtf <- sort(gtf) ## we sort all exons according to start on genome

## replace possible "." in the transcript_id names
# mcols(gtf)$transcript_id <- gsub(".", "_", mcols(gtf)$transcript_id, fixed=TRUE)

inds <- split(1:length(gtf),as.factor( mcols(gtf)$transcript_id) )  ## get exon index per transcript
inds1 <- sapply(inds,.subset,1)  ## only the first exon per transcript
gtfs <- split(gtf, as.factor( mcols(gtf)$transcript_id) ) ## list of all transcripts
ws <- width(gtfs)
m <- match(names(gtfs), names(inds1))
str <- as.character(strand(gtf))[inds1[m]] ## strand of all transcripts


getPositions <- function(u, do_rev=FALSE) {
  if (do_rev)
    cs <- cumsum(rev(u))
  else
    cs <- cumsum(u)
  z <- cbind(c(1,cs[-length(cs)]+1),cs)
  if (do_rev)
    return(z[nrow(z):1,,drop=FALSE])
  else
    return(z)
}

ps <- vector("list",length(ws))
names(ps) <- names(ws)
s <- str=="+"
ps[s] <- lapply(ws[s], getPositions)
ps[!s] <- lapply(ws[!s], getPositions, do_rev=TRUE)

mcols(gtf)$pos_start <- unsplit(lapply(ps,function(u) u[,1]), as.factor( mcols(gtf)$transcript_id) )
mcols(gtf)$pos_end <- unsplit(lapply(ps,function(u) u[,2]), as.factor( mcols(gtf)$transcript_id) )




####### alternative: sort the exons according to position in transcript (decreasing for transcripts on the "-" strand)
# gtf_plus <- gtf[strand(gtf) == "+"]
# gtf_plus <- sort(gtf_plus)
# gtf_neg <- gtf[strand(gtf) == "-"]
# gtf_neg <- sort(gtf_neg, decreasing=TRUE)
# gtf <- c(gtf_plus, gtf_neg)

# gtfs <- split(gtf, as.factor(mcols(gtf)$transcript_id))  ## list of all transcripts
# ws <- width(gtfs)

# getPositions <- function(u) {
#   cs <- cumsum(u)
#   z <- cbind(c(1,cs[-length(cs)]+1),cs)
#   return(z)
# }

# ps <- vector("list",length(ws))
# names(ps) <- names(ws)
# # s <- str=="+"
# ps <- lapply(ws, getPositions)

# mcols(gtf)$pos_start <- unsplit(lapply(ps,function(u) u[,1]), as.factor(mcols(gtf)$transcript_id) )
# mcols(gtf)$pos_end <- unsplit(lapply(ps,function(u) u[,2]), as.factor(mcols(gtf)$transcript_id))





## annotation grouped by the exon location
# quant$BoundaryLength <- quant$Length - quant$EffectiveLength

## if two times the boundary is longer than the actual transcript, divide it by two
quant$BoundaryLength <- ifelse(2*(quant$Length - quant$EffectiveLength) > quant$Length, (quant$Length - quant$EffectiveLength)/2, quant$Length - quant$EffectiveLength)


# quant$Coverage <- quant$NumReads * fragmentLength / (quant$EffectiveLength)
# quant$Slope <- quant$Coverage / quant$BoundaryLength


matched <- match(mcols(gtf)$transcript_id, quant$Name)
quant_m <- quant[matched,]
mcols(gtf)$trLength <- quant_m$Length
mcols(gtf)$effLength <- quant_m$EffectiveLength
mcols(gtf)$boundaryLength <- quant_m$BoundaryLength


# mcols(gtf)$coverage <- quant_m$NumReads * fragmentLength / (quant_m$EffectiveLength)

# mcols(gtf)$slope <-  quant_m$Coverage / quant_m$BoundaryLength
mcols(gtf)$reads <- quant_m$NumReads
mcols(gtf)$exon_loc_id <- paste0(seqnames(gtf),"_",start(gtf),":",end(gtf))


# pos_trans <- split(mcols(gtf)[,c("transcript_id", "pos_start", "pos_end")], ex) 

# gtfs <- split(gtf, mcols(gtf)$transcript_id)  ## list of all transcripts



#### get the location of the exon relative to the effectiveLength region
mid <- which(mcols(gtf)$pos_start > mcols(gtf)$boundaryLength & mcols(gtf)$pos_end <= mcols(gtf)$trLength - mcols(gtf)$boundaryLength) 
withinLeftBound <- which( mcols(gtf)$pos_start < mcols(gtf)$boundaryLength  & mcols(gtf)$pos_end <= mcols(gtf)$boundaryLength)
withinRightBound <- which( mcols(gtf)$pos_start > mcols(gtf)$trLength - mcols(gtf)$boundaryLength  & mcols(gtf)$pos_end > mcols(gtf)$trLength - mcols(gtf)$boundaryLength)
overlapLeftBound <- which( mcols(gtf)$pos_start < mcols(gtf)$boundaryLength  & mcols(gtf)$pos_end > mcols(gtf)$boundaryLength & mcols(gtf)$pos_end <= mcols(gtf)$trLength - mcols(gtf)$boundaryLength)
overlapRightBound <- which(mcols(gtf)$pos_start > mcols(gtf)$boundaryLength &  mcols(gtf)$pos_start <= mcols(gtf)$trLength - mcols(gtf)$boundaryLength & mcols(gtf)$pos_end > mcols(gtf)$trLength - mcols(gtf)$boundaryLength)
overlapBoth <- which( mcols(gtf)$pos_start <= mcols(gtf)$boundaryLength  & mcols(gtf)$pos_end > mcols(gtf)$trLength - mcols(gtf)$boundaryLength )


location <- vector("character", length=length(gtf))
location[mid] <- "mid"
location[withinLeftBound] <- "left bound"
location[withinRightBound] <- "right bound"
location[overlapLeftBound] <- "overlap left bound"
location[overlapRightBound] <- "overlap right bound"
location[overlapBoth] <- "overlap both bounds"
mcols(gtf)$location <- as.factor(location)





# ex <- paste0(seqnames(gtf),"_",start(gtf),":",end(gtf))      # exon locations pasted
# gtf_trans <- split(mcols(gtf)$transcript_id, ex)             # transcripts grouped by the exon location

# tkey1 <- with(truth, paste0(seqnames,"_",start,":",end))
# m1 <- match(names(gtf_trans), tkey1)
# truth1 <- truth[m1,] ## get the same order of exons than in gtf_trans 

## convert it to a data frame
gtf_df <- as.data.frame(mcols(gtf)[,c("transcript_id", "exon_loc_id", "pos_start", "pos_end", "trLength", "effLength", "boundaryLength", "reads", "location" )])
gtf_df$exonLength <- width(gtf)


### split exons according to their location
# gtf_exonID <- split(gtf_df, gtf_df$exon_loc_id)
# ## read the truth and match the exon from the transcript list to the exons in the truth table 
# tkey <- with(truth, paste0(seqnames,"_",start,":",end))
# m1 <- match(names(gtf_exonID), tkey1)
# truth <- truth[m1,] ## get the same order of exons than in gtf_trans 





#### compute the coverage using the possible fragment start positions   ==============================
startPos <- pmax(mcols(gtf)$pos_start - ( fragmentLength -1), 1)
endPos <- pmin(mcols(gtf)$pos_end + (fragmentLength -1), mcols(gtf)$trLength)
# nFragPos <- endPos - startPos + 1

nforwPos <- mcols(gtf)$pos_end - startPos
nrevPos <- endPos - mcols(gtf)$pos_start 
nFragPos <- pmax(nforwPos, nrevPos)

mcols(gtf)$nFragPos <- nFragPos
gtf_df$nFragPos <- nFragPos



##### write a function that takes into account the location of the exon and that does compute the scaling factor for each transcript


## The boundaryLength = transcript length - effective Length


getScalingFactor <- function(u){
    if(unique(u$location) == "mid"){ ## coverage is allready correct
        rep(1, nrow(u))
    } else if(unique(u$location) == "right bound"){
	    ## exon in right boundary
	    # compute the midpoint of the exon relative to boundary start
	    xmid <- u$pos_start + ( u$pos_end - u$pos_start) / 2
	    # position of midpoint in right boundary
	    xmid <- xmid - ( u$trLength - u$boundaryLength)
	    ## compute height at xmid  = scaling factor
	    # slope is -1/boundaryLength
	    # y = -1/F * x + 1
	    # F is boundaryLength
	    ymid <- -1/u$boundaryLength * xmid + 1	  
	    return(ymid)
	} else if(unique(u$location) == "left bound"){
	## exon is in left boundary
	    # compute the midpoint of the exon
	    xmid <- u$pos_start + ( u$pos_end - u$pos_start) / 2
	    ## compute height at xmid  = scaling factor
	    # slope is -1/boundaryLength
	    # y = 1/F * x
	    # F is boundaryLength
	    ymid <- 1/u$boundaryLength * xmid  
	    return(ymid)
	} else if(unique(u$location) == "overlap right bound"){
		## exon overlaps with the right boundary
		# compute xmid for the part that overlaps with the boundary
		xmid <- ( u$trLength - u$boundaryLength) + ( u$pos_end - ( u$trLength - u$boundaryLength)) / 2
		xmid <- xmid - ( u$trLength - u$boundaryLength)
 		ymid <- -1/u$boundaryLength * xmid + 1	

		# compute how much of the exon is located in the boundary 
		## effLen + 1 ?????
		b <- u$pos_end - ( u$trLength - u$boundaryLength) + 1
		## compute the scaling factor for the whole exon
		# sf <- (b * ymid + (u$nFragPos - b) ) / u$nFragPos
		sf <- (b * ymid + (u$exonLength - b) ) / u$exonLength
		return(sf)
	} else if(unique(u$location) == "overlap left bound"){
		## exon overlaps with the left boundary
		# compute xmid for the part that overlaps with the boundary
		xmid <- u$pos_start + (u$boundaryLength - u$pos_start) / 2
		ymid <- 1/u$boundaryLength * xmid    
		# compute how much of the exon is located in the boundary 
		b <- u$boundaryLength - u$pos_start
		## compute the scaling factor for the whole exon
		# sf <- (b * ymid + (u$nFragPos - b) ) / u$nFragPos
		sf <- (b * ymid + (u$exonLength - b) ) / u$exonLength
		return(sf)
	} else{ 
	## exon overlaps with both boundaries

		# compute xmid for the part that overlaps with the left boundary
		xmidleft <- u$pos_start + (u$boundaryLength - u$pos_start) / 2
		ymidleft <- 1/u$boundaryLength * xmidleft    
		# compute how much of the exon is located in the boundary 
		bleft <- u$boundaryLength - u$pos_start

		# compute xmid for the part that overlaps with the boundary
		xmidright <- ( u$trLength - u$boundaryLength) + ( u$pos_end - ( u$trLength - u$boundaryLength)) / 2
		xmidright <- xmidright - ( u$trLength - u$boundaryLength)
 		ymidright <- -1/u$boundaryLength * xmidright + 1	

		# compute how much of the exon is located in the boundary 
		## effLen + 1 ?????
		bright <- u$pos_end - ( u$trLength - u$boundaryLength) + 1
		
		## compute the scaling factor for the whole exon
		sf <- (bleft * ymidleft + bright * ymidright  +  (u$exonLength - bleft - bright) ) / u$exonLength
		return(sf)
	}
}


## split exons according to location
gtf_loc <- split(gtf_df, droplevels(gtf_df$location) )
## compute the scaling factor
scalingFactors <- sapply(gtf_loc, getScalingFactor)
# add it to data frame
gtf_df$scalingFactor <- unsplit(scalingFactors, droplevels(gtf_df$location))

# gtf_df$scalingFactors <- unsplit( sapply(split(gtf_df, droplevels(gtf_df$location) ), getScalingFactor), droplevels(gtf_df$location))




## compute the exon counts  ====================================
gtf_exonID <- split(gtf_df, gtf_df$exon_loc_id)

getFragPosCounts <- function(u){
    exCounts <- sum(u$reads / u$effLength *  u$nFragPos * u$scalingFactor)
    return(exCounts )
}
cov_fragPos <-  sapply(gtf_exonID, getFragPosCounts)


tkey <- with(truth, paste0(seqnames,"_",start,":",end))
pkey <- names(cov_fragPos) 

m1 <- match(pkey, tkey)
pkey <- pkey[!is.na(m1)]  ## some exons might not be annotated because they are wrong predictions
cov_fragPos <-  cov_fragPos[!is.na(m1)]
m1 <- match(pkey, tkey) 
truth <- truth[m1,] ## get the same order of exons than in gtf_trans 

type <- rep("annotated", length(cov_fragPos) )
type[pkey %in% removed_me_loc] <- "novel_me"
type[pkey %in% removed_exon_loc] <- "novel_exon"	



#### write the results table
results <- data.frame(seqnames = truth$seqnames, start = truth$start, end = truth$end, salmon_count = cov_fragPos, type = type)

write.table(results, file.path(SalmonOutput, "salmon_coverage_count.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE)








#### ------------------- FIGURES ------------------------ ###


### make individual plots for microexons and normal exons for BC2 poster

df_cov_fragPos <- data.frame(salmon_count=round(cov_fragPos, 0), 
                 exon=c("NOT microexons","microexons")[1+(truth$length<=27)], 
                 read_count=truth$count_reads, 
                 exonLength=truth$length, type = type )

### all microexons
p <- ggplot(df_cov_fragPos[df_cov_fragPos$exon=="microexons",], aes(x=salmon_count, y=read_count, color=exonLength, fill=..count..)) +
            stat_bin_hex(bins=60 ) + scale_fill_viridis()+  
            scale_x_log10() + scale_y_log10() +
            geom_abline(slope=1, intercept=0, colour="black", size=0.6) +
            theme_bw(base_size = 14) + xlab("derived exon counts") + ylab("true exon counts") + ggtitle("microexons")
ggsave(paste0(SalmonOutput, "salmon_coverageCounts_vs_readcount_nFragPos_scaled_binhex_microexons60.pdf"), p, device="pdf",  width=6, height=5,units="in")



### all normal exons
p <- ggplot(df_cov_fragPos[df_cov_fragPos$exon=="NOT microexons",], aes(x=salmon_count, y=read_count, color=exonLength, fill=..count..)) +
            stat_bin_hex(bins=60 ) + scale_fill_viridis()+  
            scale_x_log10() + scale_y_log10() +
            geom_abline(slope=1, intercept=0, colour="black", size=0.6) +
            theme_bw(base_size = 14) + xlab("derived exon counts") + ylab("true exon counts") + ggtitle("NOT microexons")
ggsave(paste0(SalmonOutput, "salmon_coverageCounts_vs_readcount_nFragPos_scaled_binhex_NOTmicroexons60.pdf"), p, device="pdf",  width=6, height=5,units="in")


## plot only the novel me 
p <- ggplot(df_cov_fragPos[df_cov_fragPos$type=="novel_me",], aes(x=salmon_count, y=read_count, color=exonLength, fill=..count..)) +
            stat_bin_hex(bins=60 ) + scale_fill_viridis()+  
            scale_x_log10() + scale_y_log10() +
            geom_abline(slope=1, intercept=0, colour="black", size=0.6) +
            theme_bw(base_size = 14) + xlab("derived exon counts") + ylab("true exon counts") + ggtitle("novel microexons")
ggsave(paste0(SalmonOutput, "salmon_coverageCounts_vs_readcount_nFragPos_scaled_binhex_60_novel_me.pdf"), p, device="pdf",  width=6, height=5,units="in")

## only the novel exons
p <- ggplot(df_cov_fragPos[df_cov_fragPos$type=="novel_exon",], aes(x=salmon_count, y=read_count, color=exonLength, fill=..count..)) +
            stat_bin_hex(bins=60 ) + scale_fill_viridis()+  
            scale_x_log10() + scale_y_log10() +
            geom_abline(slope=1, intercept=0, colour="black", size=0.6) +
            theme_bw(base_size = 14) + xlab("derived exon counts") + ylab("true exon counts") + ggtitle("novel exons")
ggsave(paste0(SalmonOutput, "salmon_coverageCounts_vs_readcount_nFragPos_scaled_binhex_60_novel_exon.pdf"), p, device="pdf",  width=6, height=5,units="in")


### all microexons, annotated and novel me are color coded
p <- ggplot(df_cov_fragPos[df_cov_fragPos$exon=="microexons",], aes(x=salmon_count, y=read_count, color = type)) +
            geom_point(alpha=.8) + 
            scale_x_log10() + scale_y_log10() +
            geom_abline(slope=1, intercept=0, colour="black", size=1) +
            theme_bw(base_size = 14) + xlab("Salmon counts") + ylab("true counts") +
            scale_color_manual(values=c("grey", "red"))
ggsave(paste0(SalmonOutput, "salmon_coverageCounts_vs_readcount_nFragPos_scaled_microexons.pdf"), p, device="pdf",  width=6, height=5,units="in")


p <- ggplot(df_cov_fragPos[df_cov_fragPos$exon=="NOT microexons",], aes(x=salmon_count, y=read_count, color = type)) +
            geom_point(alpha=.2) + 
            scale_x_log10() + scale_y_log10() +
            geom_abline(slope=1, intercept=0, colour="black", size=1) +
            theme_bw(base_size = 14) + xlab("Salmon counts") + ylab("true counts") + 
            scale_color_manual(values=c("grey", "red"))

ggsave(paste0(SalmonOutput, "salmon_coverageCounts_vs_readcount_nFragPos_scaled_exons.pdf"), p, device="pdf",  width=6, height=5,units="in")








#### TODO check why we have negative  scalingFactors!
### Todo: check that all positions are correct and we are not missing a - or + 1
#### TODO: do we have to scale for the nFragPos? or just ignore it and multiply with scalilng facgtor?
###### TODO headmap to visualize the number of points
### TODO: find correct boundary length!!!!

# head(gtf_df[gtf_df$scalingFactor < 0,])
#      transcript_id    exon_loc_id pos_start pos_end trLength effLength
# 6  ENST00000606546 19_66346:66499        59     212      236   29.1439
# 7  ENST00000606546 19_60105:60162         1      58      236   29.1439
# 16 ENST00000607633 19_66346:66382       170     206      206   22.2098
# 24 ENST00000391654 19_70928:71535       155     762      762  536.0030
# 26 ENST00000606592 19_71431:71568       155     292      292   67.6751
# 27 ENST00000606592 19_66346:66499         1     154      292   67.6751
#    boundaryLength     reads            location exonLength nFragPos
# 6        206.8561 0.0000000 overlap both bounds        154      211
# 7        206.8561 0.0000000 overlap both bounds         58      235
# 16       183.7902 0.0000000 overlap both bounds         37      205
# 24       225.9970 0.0136024 overlap both bounds        608      761
# 26       224.3249 0.0000000 overlap both bounds        138      291
# 27       224.3249 0.0000000 overlap both bounds        154      291
#    scalingFactor
# 6   -0.039012070
# 7   -0.874471111
# 16  -2.114678335
# 24  -0.090501787
# 26  -0.386603223
# 27  -0.002021906


# --> problems if transcript length ~ fragment length and if effLength is really short


##### we have 4 cases, where the exon end position is after the transcript end position!!!!
# > gtf_df[gtf_df$scalingFactor < 0,]
#         transcript_id          exon_loc_id pos_start pos_end trLength effLength
# 42504 ENST00000590828 19_36387120:36387161      1511    1552     1531   1304.09
# 44570 ENST00000589632 19_38270184:38270200      2485    2501     2486   2259.09
# 55400 ENST00000330997 19_44905734:44905755      5992    6013     5993   5766.09
# 66195 ENST00000593919 19_50528603:50528629      1953    1979     1959   1732.09
#       boundaryLength      reads    location exonLength nFragPos scalingFactor
# 42504         226.91   67.99110 right bound         42      270  -0.002203517
# 44570         226.91    7.40478 right bound         17      245  -0.030849235
# 55400         226.91 2587.97000 right bound         22      250  -0.041866819
# 66195         226.91  167.56500 right bound         27      255  -0.030849235






## maybe use switch because it is faster!

# centre <- function(x, type) {
#        switch(type,
#               mean = mean(x),
#               median = median(x),
#               trimmed = mean(x, trim = .1))
#      }
#      x <- rcauchy(10)
#      centre(x, "mean")
#      centre(x, "median")
#      centre(x, "trimmed")


# ccc <- c("b","QQ","a","A","bb")
#      # note: cat() produces no output for NULL
#      for(ch in ccc)
#          cat(ch,":", switch(EXPR = ch, a = 1, b = 2:3), "\n")
#      for(ch in ccc)
#          cat(ch,":", switch(EXPR = ch, a =, A = 1, b = 2:3, "Otherwise: last"),"\n")

