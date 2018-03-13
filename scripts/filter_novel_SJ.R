#### This script takes the STAR output (SJ.tab.out) and looks for 2 consecutive SJs within an annotated SJ.
## File rlist based on # of reads
##						# of sampels in which the SJ occur (at least x reads in all samples, ...)
## --> output: list of novel splicing events with start + stop coordinates and end of previous and start of following exon

library(rtracklayer)
library(data.table)
library(GenomicFeatures)
library(GenomicAlignments)

library(dplyr)


GTF <- snakemake@input[["gtf"]]
SJFILE <- snakemake@input[["sj"]]
BAM <- snakemake@input[["bam"]]
OUTFILE <- snakemake@output[["outfile"]]


# ## Defining file paths
# GTF="simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf"
# SJFILE="simulation/mapping/STAR/me/outSJfilterOverhangMin6/pass2_SJ.out.tab"
# OUTFILE="test_novel_exons/me_outSJfilterOverhangMin6.txt"
# BAM <- "simulation/mapping/STAR/me/outSJfilterOverhangMin6/pass2_Aligned.out_s.bam"

# GTF <- "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf"
# SJFILE <- "simulation/mapping/STAR/me/outSJfilterOverhangMin9/pass2_SJ.out.tab"
# BAM <- "simulation/mapping/STAR/me/outSJfilterOverhangMin9/pass2_Aligned.out_s.bam"
# OUTFILE <- "simulation/analysis/filtered_SJ/novel_exons_reduced_me_outSJfilterOverhangMin9.txt"

# GTF <- "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_exon.gtf"
# SJFILE <- "simulation/mapping/STAR/exon/outSJfilterCountTotalMin3/pass2_SJ.out.tab"
# BAM <- "simulation/mapping/STAR/exon/outSJfilterCountTotalMin3/pass2_Aligned.out_s.bam"
# OUTFILE <- "simulation/analysis/filtered_SJ/novel_exons_reduced_exon_outSJfilterCountTotalMin3.txt"

# GTF <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me_exon.gtf"
# SJFILE <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/mapping/STAR/me_exon/outSJfilterOverhangMin6/pass2_SJ.out.tab"
# BAM <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/mapping/STAR/me_exon/outSJfilterOverhangMin6/pass2_Aligned.out_s.bam"



## Reading files
print("Reading data")

gtf <- import(GTF)
sj <- fread(SJFILE)
# colnames(sj) <- c("chr", "first", "last", "strand", "motif", "annotated", "unique", "mutimapping", "maxoverhang")
colnames(sj) <- c("seqnames", "start", "end", "strand", "motif", "annotated", "unique", "mutimapping", "maxoverhang")

#strand: (0: undefined, 1: +, 2: -)
sj$strand[sj$strand==0] <- "*"
sj$strand[sj$strand==1] <- "+"
sj$strand[sj$strand==2] <- "-"

##########################
# junction filtering
#########################


# filter out all sj with annotated==0, because the 2nd pass mapping is only for quantification and not for discovery (they are most likely artifacts!
## 0: unannotated, 1: annotated 
sj <- sj[sj$annotated != 0, ]  ## remove 6 junctions
## but this filters out one of the true junctions that we want to find!!

## Extract introns from the gtf file

## load gene model from GTF file
exons <- gtf[mcols(gtf)$type =="exon", ]
txdb <- makeTxDbFromGRanges(exons)
# txdb <- makeTxDbFromGFF(GTF, format="gtf")
## get all annotated introns 
inbytx <- intronsByTranscript(txdb, use.names=TRUE)
introns <- unique(unlist(inbytx))

## overlap introns and SJ to filter out all annotated SJs
## we cannot use the "annotated" column, because in the 2-pass alignment, all SJ from the first pass are treated as annotated, even if they are not in the GTF file!
print("Filter SJ")

# sj.gr <- GRanges(seqnames=sj$chr, ranges=IRanges(start=sj$first, end=sj$last), strand=sj$strand, motif=sj$motif, annotated=sj$annotated, unique=sj$unique, mutimapping=sj$mutimapping, maxoverhang=sj$maxoverhang)

sj.gr <- GRanges(sj)

## filter out all annotated junctions
sj.ann <- subsetByOverlaps(sj.gr, introns, type="equal") 
sj.unann <- sj.gr[!(sj.gr %in% sj.ann)]

## filter out all junctions with less than XXXX reads
# sj.unann <- mcols(sj.unann)$unique >= 1


## filter out all novel junctions that do not touch an annotated exon on the same strand
print("Filter sj on wrong strand")
sj.unann <- sj.unann[(start(sj.unann)-1) %in% end(exons) | (end(sj.unann)+1) %in% start(exons),]

## Which side of the sj is touched by an exon?
## Find all exons with the same strand as the splice junction and that touch it
## if both sj touch an exon, the exons must come from the same gene
## return: a string that defines which side of the sj is toughing an exon: both sides, start or end of sj
sj_touching_exon <- function(sj, exons){
  s <- which((start(sj)-1) == end(exons) ) 
  s <- exons[s,]   ## all touching exons at start or sj
  s <- s[strand(s) == strand(sj),] ## all touching exons with same strand

  e <- which((end(sj)+1) == start(exons))  
  e <- exons[e]   ## all touching exons at end of sj
  e <- e[strand(e) == strand(sj),]

if( length(s) >0 ){
  if(length(e) >0){
    if(any( mcols(s)$gene_id %in% mcols(e)$gene_id )) { ## pair of touching exons with same strand
    "both"
    } else NA
  } else "start"
} else if(length(e) >0){
  "end"
  } else NA

}

touching <- sapply(sj.unann, function(x) sj_touching_exon(x, exons=exons))
mcols(sj.unann)$touching <- touching ## either start or end, NA if not known
sj.unann <- sj.unann[!is.na(touching)] 

### is the sj internal or terminal?
mcols(sj.unann)$location <- NA ## internal or terminal, Na if not known

## Look for two consecutive novel splice junctions that are located within an annotated junction, but have the same start or end --> the 2 novel SJ could replace the annotated intron 
##  XXX--------XXX     annotation
##  XXX---X----XXX     novel found exon 

## First we get all SJ that are located within annotated junctions
## Then, we select only the junctions that have the same start or end as the intron
## Last, we need to find all introns that have both a novel start and end SJ.
## The novel SJs have to be consecutive.

print("Find sj pairs in same intron")

# within <- subsetByOverlaps(sj.unann, introns, type="within")
within_ind <- unique(queryHits(findOverlaps(sj.unann, introns, type="within")))
within <- sj.unann[within_ind,]

startHit <- findOverlaps(within, introns, type="start")
endHit <- findOverlaps(within, introns, type="end")
mcols(sj.unann)$location[within_ind[ c(queryHits(startHit), queryHits(endHit))]] <- "internal"  ## all sj that connect with their start or end to an annotated intron are "internal" splice junctions 


## Try to pair sj in the same intron
startHit <- startHit[subjectHits(startHit) %in% subjectHits(endHit), ] ## filter out introns with only start/end
endHit <- endHit[subjectHits(endHit) %in% subjectHits(startHit), ]

starts <- split(within[queryHits(startHit)], as.factor(subjectHits(startHit)))  ## get GRanges per intron
ends <- split(within[queryHits(endHit)], as.factor(subjectHits(endHit)))

## This function matches the novel SJ within one intron and makes sure that they are consecutive
filterSJ <- function(s, e){
  if(all(length(s) >1, length(e)>1)){
    print("More than one novel SJ.")
  }
  else{
    if(isDisjoint(c(s, e), ignore.strand=FALSE)){
      return(c(as.vector(seqnames(s)), start(s)-1, end(s)+1, start(e)-1, end(e)+1, as.vector(strand(s)) ) )
    }else{
      return(NA)
    }
  }
}

if(length(starts) > 0 & length(ends) > 0 ){
  novelExons <- lapply(names(starts), function(x) filterSJ(s = starts[[x]], e = ends[[x]]))


  ### the introns that returned NA are the introns where the two splice junctions were not consective
  ### --> we keep these junctions, they might be terminal ones or cases where one end ot the exon is already annotated
  matched_starts <- unlist(starts[!is.na(novelExons)]) ## the internal start junctions
  matched_ends <- unlist(ends[!is.na(novelExons)]) ## the internal start junctions
  sj.unann <- subsetByOverlaps(sj.unann, c(matched_starts, matched_ends), type="equal", invert=TRUE)

  novelExons <- as.data.frame(matrix(unlist(novelExons), ncol=6, byrow=T), stringsAsFactors=FALSE )
  colnames(novelExons) <- c("seqnames", "lend", "start", "end", "rstart", "strand")
} else{
  novelExons <- data.frame(seqnames= character(), lend = integer(), start = integer(), end = integer(), rstart = integer(), strand = character(), stringsAsFactors=FALSE )
}

###### TODO: why is the second cassette exon of 19 56180810 56181058 not found???
### is 56180547 56180810 not a novel SJ???

#    seqnames     lend    start      end   rstart strand
# 25       19 56180535 56180810 56181058 56185300      +




### find possible terminal exons
### read in all reads that overlap with the splice junctions and find out if it is a terminal or internal exon
### i.e. are there reads with the splice junction, but also a second splice junctions that tells us the size of the novel exon?

## read all reads that have mapped nucleotides in the junction region (including last and first nucleotide of the two exons)
print("Read in BAM file and filter all reads with the novel sj")
## read in all reads
param <- ScanBamParam(which = sj.unann, what = c("qname"))
reads <- readGAlignments(BAM, index = BAM, with.which_label=TRUE, param=param) ## all reads that overlap 
## filter out all reads with 0 junctions  --> this is still based on the pairs!!
reads <- reads[njunc(reads) > 0,]
## filter out reads that do not contain the junction with which they overlap
junc <- junctions(reads, use.mcols=TRUE) ## list with all junctions from each read
# names(junc) <- mcols(junc)$which_label
names(junc) <- 1:length(junc)
true_junc <- as.character(mcols(junc)$which_label) ##the unannotated junctions

junc <- unlist(junc)

df <- data.frame(read_nr = as.integer(names(junc)), junction = paste0(seqnames(junc), ":", start(junc), "-", end(junc)), true_junction = rep(true_junc, njunc(reads)), ## number of junctions per reads
  stringsAsFactors=FALSE)

reads <- reads[ df$read_nr[df$junction == df$true_junction] ] ## keep all reads that contain the novel junction
rm(df)



#### find exon coordinates for all remaining sj: either terminal or overlapping annotated exons

#### get the start and end position of the transcript
## This is only needed for the simulated data, because the transcript annotation still contains the removed exons in the GTF file. For real data, we can simply take to coordinates of the entry with type="transcript".
transcript_range <- function(gr){
  GRanges(seqnames = seqnames(gr)[1], ranges = IRanges(min(start(gr)), max(end(gr))), strand= strand(gr)[1])
}

#######This function computes if and and which exon of a sj is terminal:
## overlap sj with all genes
## overlap with all transcripts of the gene:
## take all exons per transcript and comput the start and end positions
## remove all overlappping transcripts
## take all exons of the remaining transcripts
## find the exons where the sj spliced to
## either "start" or "end" of the sj

## j: novel SJ
## txdb: txdb object from the GTF file
## gtxdb: all genes from the txdb
## ebyTr: all exons per transript
which_exon_terminal <- function(j, txdb, gtxdb, ebyTr){
  ### maybe we have to include the position +-1 of the sj, in case of terminal sj that is outside of the annotated gene boundary ??
  genes <- mcols( subsetByOverlaps(gtxdb, j) )$gene_id ## all genes that overlap with the SJ
  tr <- transcripts(txdb, filter = list(gene_id = genes))
  tr <- mcols( tr )$tx_name  ## all transcripts from the genes
  e <- ebyTr[tr]
  tr_range <- unlist(GRangesList(lapply(e, transcript_range)))  ## the range of all transcripts
  tr <- names( subsetByOverlaps(tr_range, j, invert = TRUE) ) ## all transcripts that do NOT overlap with the sj
  # tr <- subsetByOverlaps(transcripts(txdb, filter = list(gene_id = genes)), j, invert = TRUE) ## Use this for real data, where the transcript boundaries correspond to the start and end of the first and last exon, but it does not work for the simulated data, because the "transcript" entry in the gtf file still contains the remvoved exons!
 
  e <- unlist(ebyTr[tr]) ## all exons from the non-overlapping transcripts
  if(length(e) == 0){print("Sj overlaps with all transcripts of the gene."); return(NA)}

  ## does the sj touch any of the exons?
  start <- start(j)-1 == end(e)
  end <- end(j) + 1 == start(e)

  return(ifelse(any(start), "start", ifelse(any(end), "end", "NA")))
}


## This function returns the end of the novel exon and the start of the consecutive exon
## or only the end if the exon is terminal
## if touching == "start"
## j: junction 
## x: mapped ranges of read
get_exon_end_coordinate <- function(j, x){
  exon_id <- which( (end(j)+1) == start(x) )
  if(exon_id < length(x)) c(end(x[exon_id]), start(x[exon_id+1]) ) else c(end(x[exon_id]), NA)
}

## This function return the end of the previous exon and the start of the novel exon
## or only the start if the exon is terminal
## if touching =="end"
## j: junction 
## x: mapped ranges of read
get_exon_start_coordinate <- function(j, x){
   exon_id <- which( (start(j)-1) == end(x) )
   if(exon_id>1) c(end(x[exon_id-1]), start(x[exon_id])) else c(NA, start(x[exon_id]) )
}



## Determine the size of the novel exon (un-paired internal sj and terminal junctions) based on the touching exons and the reads:
## touching == start or end:
##    check if there are any reads with two sj --> we know the size of the exon and where it splices to
## touching == both:
##    complicated overlapping exons: either terminal, or overlapping exon annotations
##    --> still check the reads, but if they do not help because it is terminal, return NA
## 
get_exon_boundary <- function(j, reads){
  j_string <- toString(j)
  j_string <- substring(j_string, 1, nchar(j_string)-2) # without the strand, because it is missing in the which_label of the reads and because the read strand is independent from sj strand ("+" forward and "-" reverse read)
  r <- reads[ mcols(reads)$which_label == j_string ] ## all reads with the sj

  r_mapped <- grglist(r, order.as.in.query=FALSE) ##this are the mapped ranges of each read, from left to right (independent of strand)
  # r_mapped <- r_mapped[ unlist(lapply(r_mapped,length)) > 2 ]## all reads with at least 2 sj

  # rj <- junctions(r)
  #   ## filter out all reads with only one junction
  # rj_multi <- unlist(lapply(rj,length)) > 1
  # rj <- rj[rj_multi]

  
  if(mcols(j)$touching == "start"){  ## the start of the novel SJ touches an annotated exon
  ##    X---X        novel junction
  ##  xxx---xx----x  read
  ##         X----X  want to find 
  coord <- unique(t(sapply(r_mapped, function(x) get_exon_end_coordinate(j = j, x = x)) ))
  ## if there are reads with a second splice junctions, we know the start coordinates of the consecutive exon
  ## if not, then we do not even know the end of the novel exon, so we simply take the lngest read --> max end coordinatecd
  c(as.vector( seqnames(j) ), 
    if( any( !is.na(coord[,2]) )) c(start(j)-1, end(j)+1, coord[!is.na(coord[,2]),] ) 
    else c(start(j)-1, end(j)+1, max(coord[,1]), NA), 
    as.vector(strand(j)) )
  }else if(mcols(j)$touching == "end"){
    ##         X----X   novel junction
    ##   xx---xx----xx  read
    ##    X---X         want to find
    coord <- unique( t(sapply(r_mapped, function(x) get_exon_start_coordinate(j = j, x = x)) ) )
    ## matrix with two columns: lend start
    ## if there are reads with a second splice junctions, we know the end coordinates of the previous exon
    ## if not, then we do not even know the start of the novel exon, so we simply take the longest read --> min start coordinate
    c(as.vector( seqnames(j) ), 
      if( any( !is.na(coord[,1]) ) ) c(coord[!is.na(coord[,1]),], start(j)-1, end(j)+1 ) 
      else c(NA, min(coord[,2]), start(j)-1, end(j)+1 ), 
      as.vector( strand(j) ) )
  }else{
    ##     X---X       novel junction
    ##         xxxxx   annotated exon
    ##  xxxx           annotated exon    
    ##    nn---xxxxx   possible transcript with novel exon nn at the start of the SJ
    ##  xxxx---nn      possible transcript with novel exon nn at the end of the SJ

    ## We do not know if there is a novel exon and if yes, which side of the novel sj it connects to.
    ## The exon is most likely a terminal exon, because otherwise we would have identified it as an internal exon
    ## --> look at the reads, maybe they have two splice junctions and help us to determine the boundaries of the novel exon.
     
    ## Try to find the coordinates of the exon at the start of the sj:
    start_coord <- unique( t(sapply(r_mapped, function(x) get_exon_start_coordinate(j = j, x = x)) ) )
    # start_coord <- start_coord[ !is.na(start_coord[,1]) , ,drop=FALSE]

    ## Try to find coordinates of the exon at the end of the sj
    end_coord <- unique(t(sapply(r_mapped, function(x) get_exon_end_coordinate(j = j, x = x)) ))
    # end_coord <- end_coord[ !is.na(end_coord[,2]),  ,drop=FALSE]

    ## If any of the two returned two coordinates, we take these as the exon prediction
    if(any( !is.na(start_coord[,1]) ) ){
      if(any( !is.na(end_coord[,2]) )){ ## we found exon coordinates at both ends of the sj
        ### TODO only return the exon combinations that are not annotated in a transcript yet
      rbind(   
        c(as.vector( seqnames(j) ), 
          c(start_coord[ !is.na(start_coord[,1]) ,], start(j)-1, end(j)+1 ), 
          as.vector( strand(j) ) ),
        c(as.vector( seqnames(j) ), 
          c(start(j)-1, end(j)+1, end_coord[ !is.na(end_coord[,2]), ] ),
          as.vector( strand(j) ) )
        )
      } else { ## novel exon at the start of the SJ
      c(as.vector( seqnames(j) ), 
          c(start_coord[ !is.na(start_coord[,1]) ,], start(j)-1, end(j)+1 ), 
          as.vector( strand(j) ) )
      }
    } else if(any( !is.na(end_coord[,2]) )){  ## novel exon at the end of the SJ
       c(as.vector( seqnames(j) ), 
          c(start(j)-1, end(j)+1, end_coord[ !is.na(end_coord[,2]), ] ),
          as.vector( strand(j) ) )
    } else { ## We did not find a novel exon yet
      ## Check if it could be a terminal exon --> is any of the touching exons the first or last in the transcript? 
      ## --> if yes, the SJ is terminal -->  determine the exon end from the max read end
      ## if not terminal, then return NA, because there is no way that we can determine the coordinates of the touching exons 
      ## --> error in reduction of GTF?

      ## is on of the exons terminal?
      terminal <- which_exon_terminal(j, txdb = txdb, gtxdb = gtxdb, ebyTr=ebyTr)

      if(is.na(terminal)){ ## We do not have enough information to infer the exon coordinates, ignore this novel SJ
        c(as.vector( seqnames(j) ), 
          c(NA, NA, NA, NA), 
          as.vector( strand(j) ) )
      } else if(terminal == "start"){
        ## The start of the SJ is connected to the transcript --> terminal exon is at the end of the sj
        ## take maximal exon size from end_coord
        # start_coord <- start_coord[ !is.na(start_coord[,1]) , ,drop=FALSE]
        c(as.vector( seqnames(j) ), 
          c(start_coord[ which.min(start_coord[,2]), ], start(j)-1, end(j)+1 ), 
          as.vector( strand(j) ) )
      } else if(terminal == "end"){
        ## The end of the sj is connected to transcript --> terminal exon is at the start of the sj
        ## take maximial exon size from start_coord
        # end_coord <- end_coord[ !is.na(end_coord[,2]),  ,drop=FALSE]
        c(as.vector( seqnames(j) ), 
          c(start_coord[ which.min(start_coord[,2]), ], start(j)-1, end(j)+1 ), 
          as.vector( strand(j) ) )
      } else{
        c(as.vector( seqnames(j) ), 
          c(NA, NA, NA, NA), 
          as.vector( strand(j) ) )
      }
    }
  }
}

print("Find exon coordinates of exons adjacent to novel sj")

ebyTr <- exonsBy(txdb, by = "tx", use.names = TRUE)  ## exons per transcript
gtxdb <- genes(txdb)  ## all gene ranges

######## determine the exon boundaries of all remaining sj (terminal or novel exon overlapping annotated one)
res <- as.data.frame( t( sapply(sj.unann, function(x) get_exon_boundary(x, reads= reads)) ), stringsAsFactors=FALSE   )
colnames(res) <- c("seqnames", "lend", "start", "end", "rstart", "strand")
## filter out all cases where we only have NA
res <- res[ rowSums(is.na(res[,c("lend", "start", "end", "rstart")])) <4, ]

if(nrow(res)>0){
  ## combine the results from the internal sj with this
  novelExons <- rbind( novelExons, res)
} 



########## ----------------- Testing ----------------------------------
### TODO: make this optional so that the user can decide if he wants to have better predictions, but slower performance
# samtools view -h simulation/mapping/STAR/me_exon/outSJfilterOverhangMin6/pass2_Aligned.out_s.bam | awk '{if ($6 ~ /N/) print $1, $3, $4, $6}'

### Use reads with 2 junctions or raid pairs with each one junction to predict novel exons with novel combinations of splice junctions (not annotated)
## We filter all reads with "N" in the CIGAR string
## and that are properly paired (-f 2)
junc_reads <- fread(paste0("samtools view  -f 2 ", BAM, " | awk '{if ($6 ~ /N/) print $1, $2, $3, $4, $6}'"))
names(junc_reads) <- c("qname", "flag", "seqnames", "pos", "cigar")

## We only keep all reas that are properly paired 
# junc_reads <- junc_reads[ bamFlagTest(junc_reads$flag, "isPaired") & bamFlagTest(junc_reads$flag, "isProperPair") ,]

## We infer the read strand from the flag: 
# our data comes from Illumina HiSeq 2000 (stranded TruSeq libary preparation with dUTPs )
## Thus, the the last read in pair determines the strand of the junction --> reverse the strand of the first reads
junc_reads <- cbind(junc_reads, bamFlagAsBitMatrix(junc_reads$flag, bitnames=c("isMinusStrand", "isFirstMateRead")) )
## isMinusStrand==0 & isFirstMateRead==0 --> "+"
## isMinusStrand==0 & isFirstMateRead==1 --> "-"
## isMinusStrand==1 & isFirstMateRead==0 --> "-"
## isMinusStrand==1 & isFirstMateRead==1 --> "+"
junc_reads[,strand:=ifelse(isMinusStrand + isFirstMateRead==0, "+", 
    ifelse(isMinusStrand + isFirstMateRead ==1, "-", "+"))]


## 1)  keep all reads with 2 junctions and also look for novel junction combinations
## filter all reads with 2 "N" in cigar
library(stringr)
two_junc_reads<- junc_reads[str_count(junc_reads$cigar, "N")==2,]

## we remove duplicate reads with the same junctions
two_junc_reads <- two_junc_reads[!duplicated(two_junc_reads[,-"qname", with=FALSE]),]

cigar_junc <- cigarRangesAlongReferenceSpace(cigar = two_junc_reads$cigar, flag = two_junc_reads$flag, pos = two_junc_reads$pos, ops = "N")
names(cigar_junc) <- 1:nrow(two_junc_reads) ## index of the read
# two_junc_reads$qname

cigar_junc<- as.data.table( unlist(cigar_junc))
cigar_junc$names <- as.numeric(cigar_junc$names)
# m <- match(cigar_junc$names, two_junc_reads$qname)
cigar_junc[, seqnames := two_junc_reads$seqnames[cigar_junc$names]]
cigar_junc[, strand := two_junc_reads$strand[cigar_junc$names]]
cigar_junc[, seqnames:= as.integer(seqnames)]

## the "names" indicates the original read --> we want to find novel combinations of introns within a read

## we join the read junctions with the annotated introns 
## start is the first nucleotid in the intron and end the last nucleotid (excluding exons)
# inbytx <- intronsByTranscript(txdb, use.names=TRUE)
## unique list of all transcripts and introns
intr_gtf <- unlist(inbytx)
intr_gtf$transcript_id <- names(intr_gtf)
intr_idx <- paste0(start(intr_gtf),":", end(intr_gtf),"_", mcols(intr_gtf)$transcript_id)
intr_gtf<- intr_gtf[match(unique(intr_idx), intr_idx)] ##only keep the unnique junction-transcript_id combinations

## we compare the read junctions to annotated transcript introns
# olap <- findOverlaps(cigar_junc_gr, intr_gtf) ## we can run it on the list?
olap <- findOverlaps(GRanges(cigar_junc), intr_gtf, type ="equal") 

## test is TRUE if there is a transcript that contains both junctions
# tmp = data.frame(read = queryHits(olap), tr = names(intr_gtf[subjectHits(olap)])) %>% group_by(read,tr) %>% count(.) %>% group_by(read) %>% summarise(test=any(n==2))

tmp <- data.frame(read = cigar_junc$names[queryHits(olap)], tr = names(intr_gtf[subjectHits(olap)])) %>% group_by(read,tr) %>% count(.) %>% group_by(read) %>% summarise(test=any(n==2))

tmp %>% filter(test==F) ## all reads that have unannotated junction combinations: 2,297
tmp %>% filter(test==T)  ## all reads that have annotated junction combinations: 4,964

read_id <- tmp %>% filter(test==F) %>% pull(read)

novel_reads <- cigar_junc[names %in% read_id]

##### predict the novel exons based on the novel junctions
novel_reads <- novel_reads[order(novel_reads$start),] ## sort according to junction start
read_pred <- novel_reads %>% group_by(names) %>% dplyr::summarise(lend = nth(start-1, 1), start1 = nth(end+1, 1), end1 = nth(start-1, 2), rstart = nth(end+1, 2), seqnames = unique(seqnames), strand = unique(strand))
read_pred <- read_pred %>% select(-names) %>% rename(start = start1, end = end1) %>% unique()



## convert the columns to integer instead of character
id <- c("seqnames", "lend", "start", "end", "rstart")
novelExons[,id] <- as.integer(unlist(novelExons[,id]))


### add the novel exons to the predictions from SJ
novelExons <- full_join(novelExons, read_pred)  ## predictions from both the reads and the SJ.out.tab



#### lapply --> slow
## check if the read junctions are annotated and if yes, if they occur in a single transcript
## TRUE: the read contains a novel exon, where the junctions are annoated
## FALSE: the jucntiosn are not annotated or the exon is already known
# novel_junc_comb <- function(junc, intr){
#   olap <- findOverlaps(junc, intr)
#   if(length(olap)>0){ 
#   # tr <- names(intr[subjectHits(findOverlaps(junc, intr))])
#   # !any(duplicated(tr))  ## TRUE: the junction combination is novel; FALSE: the junction combination is annotated
#  !any(duplicated(names(intr[subjectHits(findOverlaps(junc, intr))])))
#   }else{ ## the junctions are not annotated  TODO: do we want to predict these as novel exons??
#     FALSE
#   }

# }

# # junc <- cigar_junc_gr[[26]]
# # intr <- intr_gtf
# cigar_junc_gr <- split(GRanges(cigar_junc), cigar_junc$names)
# lapply(cigar_junc_gr[1:100], function(x) novel_junc_comb(x, intr_gtf))



####### using table joins:
# introns_tab <- as.data.table(intr_gtf)
# introns_tab[, seqnames:= as.integer(as.character(seqnames))]
# introns_tab[,transcript_id:=names(intr_gtf)]

# setkeyv(cigar_junc, c("start", "end", "seqnames", "strand", "width"))
# setkeyv(introns_tab, c("start", "end", "seqnames", "strand", "width"))

# merge(cigar_junc,introns_tab, all.x=TRUE, allow.cartesian=TRUE) ## left join
# merge(cigar_junc, introns_tab, nomatch=0)  ## inner join, only keep all junctions that are annotated
# test_merge<- merge(test, introns_tab, all = TRUE)## full join



## 2) keep all read pairs with each one junction --> look for novel junction combinations







##########------------ end Testing ----------------------------------

### Add columns with the number of reads supporting the left and right splice junction and the minimum of both 

## add the read counts for the left and right junction
novelExons <- novelExons %>% left_join(select(sj, seqnames, start, end, strand, unique) %>% mutate(start = start-1, end = end+1), by = c("seqnames", "lend" = "start", "start" = "end", "strand") ) %>% rename(unique_left =unique  )
novelExons <- novelExons %>% left_join(select(sj, seqnames, start, end, strand, unique) %>% mutate(start = start-1, end = end+1), by = c("seqnames", "end" = "start", "rstart" = "end", "strand") ) %>% rename(unique_right =unique  )

## take the minimum of both
novelExons$min_reads <- pmin(novelExons$unique_left, novelExons$unique_right, na.rm = TRUE)  


## Save the novel found junctions to file (start1, end1, start2, end2)
# write.table(novelExons, OUTFILE, row.names=FALSE, quote=FALSE, sep="\t")

print("Write out results")
zz<-file(description=OUTFILE,"w") 
write.table(novelExons, zz, row.names=FALSE, quote=FALSE, sep="\t"  ) 
close(zz) 










## the true transcript is ENST00000589006
# 19 [5657336, 5664147]      +   ## the transcript starts too early because I did not change it when I removed exons
## novel splic junction
# 19 [5657359, 5661498]      +

### Determine the length of the terminal exon

## special case that looks lika terminal exon, but is not, because the second junction is already annotated:
## for cases where the me overlaps with another exon with the same start, we take all reads with the sj, and if the reads have a second sj, we know that me is actually not terminal, but just overlapping an annotated sj!!
## XX---X----XXXX  right sj is novel, left one is already annotated
## XX--------XXXX
## XX---XXXXXXXXX
## this case is actually really common!!! (sj.unann 2 and 3)
# 19   [2097955, 2098050]      +
# 19   [5262984, 5265018]      -




## TODO!!!
#### There are some splice junctions that splice on both ends to annotated exons!
### are this cases where the exons are annotated, but not the transcripts?
## e.g. sj.unann[4]
# 19 [5657359, 5661498]      +



# ## exons ending in sj start
# exons[(start(sj.unann[4])-1) == end(exons)]
# 19 [5657252, 5657358]      +
# 19 [5657252, 5657358]      +
# ## ENST00000589863   ENST00000292123   ENST00000592224  ENST00000454510

# ## exons startin in sj end
# exons[(end(sj.unann[4])+1) == start(exons)]
# 19 [5661499, 5661819]      + #ENST00000589006  exon 2

# # the exons in the two sets do not share the same transcipt!!


# subsetByOverlaps(sj.unann[4], introns, type="equal") 
# findOverlap(start(sj.unann[4])-1, end(exons)



### cconvert the cigar strings to ranges
# cigarRangesAlongReferenceSpace()
# extractAlignmentRangesOnReference

# countCompatibleOverlaps() 
# grglist(reads) ## list of granges per read pair

# junc_summary <- summarizeJunctions(gal)  ## number of reads that contain each junction


