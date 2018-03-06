### This script adds the predicted novel exons and adds them to the reduced gtf file.
library(rtracklayer)
library(data.table)


GTF <- snakemake@input[["gtf"]]
EXONPRED <- snakemake@input[["exon_prediction"]]
OUTFILE <- snakemake@output[["outfile"]]


# GTF <- "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me_exon.gtf"
# GTF <- "simulation/reduced_GTF/GRCh37.85_chr19_22_reduced_me.gtf"
# EXONPRED <- "simulation/analysis/filtered_SJ/me_exon/novel_exons_outSJfilterCountTotalMin3.txt"
# EXONPRED <- "simulation/analysis/filtered_SJ/me/novel_exons_outSJfilterCountTotalMin3.txt"
# EXONPRED <- "simulation/analysis/filtered_SJ/me_exon/novel_exons_outSJfilterOverhangMin6.txt"



# OUTFILE <- "simulation/reduced_GTF_with_predicted_exons/GRCh37.85_chr19_22_reduced_me_exon_Novel_exons_reduced_me_exon_outSJfilterCountTotalMin3_version2.gtf"
# OUTFILE <- "/home/Shared/kathi/microexon_pipeline/simulation/reduced_GTF_with_predicted_exons/GRCh37.85_chr19_22_reduced_me_exon_Novel_exons_reduced_me_exon_outSJfilterCountTotalMin3.gtf"

gtf <- import(GTF)
exons <- gtf[ mcols(gtf)$type == "exon", ]
exon_pred <- fread(EXONPRED)


## Find all transcripts with exon_boundaries that splice to the novel exon
get_transcripts <- function(pred, exons){
  terminal <- 0
  ## if the lend is NA, this means the exon is terminal
  if(!is.na(pred[2]) ){
  ## find all exons that end in lend
    id1 <- mcols( subsetByOverlaps(exons, GRanges(seqnames = pred[1], ranges = IRanges(as.numeric(pred[2]), as.numeric(pred[2]) ), strand = pred[6]), type = "end", ignore.strand = FALSE) )$transcript_id
  } else{
    terminal <- 1
  }

  if(!is.na(pred[5])){
    ## find all exons that start in rstart
    id2 <- mcols( subsetByOverlaps(exons, GRanges(seqnames = pred[1], ranges = IRanges(as.numeric(pred[5]), as.numeric(pred[5])), strand = pred[6]), type = "start", ignore.strand = FALSE) )$transcript_id
  } else{
    terminal <- 2
  }


  if(terminal == 0){ ## it is a cassette exon
    ## make sure that the transcripts have no overlapping exon in the region where the new one will be located
    tr <- intersect(id1, id2)
    if(length(tr) == 0 ) NA

    olap_tr <- mcols( subsetByOverlaps( exons[ mcols(exons)$transcript_id %in% tr ], GRanges(seqnames = pred[1], ranges = IRanges(as.numeric(pred[2]) + 1, as.numeric(pred[5])-1  ), strand = pred[6]), ignore.strand = FALSE) )$ transcript_id
    return( tr[! tr %in% olap_tr] ) ## all transcripts that do not have an exon in the region where the new exon will be located
  } else if( terminal == 1) {  ## 5' terminal exon 
  ###### TODO: make sure that we only keep transcripts that start with the overlapping exons!!

    tr <- unique(id2)
    if(length(tr) == 0 ) NA

    olap_tr <- mcols( subsetByOverlaps( exons[ mcols(exons)$transcript_id %in% tr ], GRanges(seqnames = pred[1], ranges = IRanges(as.numeric(pred[3]), as.numeric(pred[5])-1 ), strand = pred[6]), ignore.strand = FALSE) )$ transcript_id   
    return(tr[! tr %in% olap_tr] )
  } else{  ## 3' terminal exon

    tr <- unique(id1)
    if(length(tr) == 0 ) NA
    ## remove all transcripts that have an exon in the connected intron or the novel exon
    olap_tr <- mcols( subsetByOverlaps( exons[ mcols(exons)$transcript_id %in% tr ], GRanges(seqnames = pred[1], ranges = IRanges(as.numeric(pred[2])+1, as.numeric(pred[4]) ), strand = pred[6]), ignore.strand = FALSE) )$ transcript_id   
    return(tr[! tr %in% olap_tr] )
  }

 }


tr_ids <- apply(as.data.frame(exon_pred) , 1, function(x) get_transcripts(x, exons))  ## transcript IDs for each novel exon



##### INSERT
##### TODO: how to fix the terminal cases???
# pred <- unlist( exon_pred[25,] )

# ### 5' terminal exon: exon_pred[51,]  gets inserted in intron but it is a terminal exon!!!
# pred <- unlist( exon_pred[51,] )
# olap_tr <- mcols( 
#   subsetByOverlaps( exons[ mcols(exons)$transcript_id %in% tr ], GRanges(seqnames = pred[1], ranges = IRanges(as.numeric(pred[3]), as.numeric(pred[5])-1 ), strand = pred[6]), ignore.strand = FALSE) 
#   )$ transcript_id   

### remove all transcripts with exon starts
##### INSERT END



## Make a new transcript with the novel exon
## copy all entries (only exons?)
## take one exon and replace the location
## change the transcript_id to ID_me_start:end
new_transcript <- function(pred, tr_ids, exons){
  type <- ifelse( (pred$end-pred$start+1) <28, "me", "exon")
  new_entries <- GRanges()
  
  tr_copy <- exons[ which( mcols(exons)$transcript_id  %in% tr_ids)]
  # mcols(tr_copy)$"exon_id" <- paste0( mcols(tr_copy)$"exon_id", "_", type, "_", pred$start, ":", pred$end )
  mcols(tr_copy)$"transcript_id" <- paste0( mcols(tr_copy)$"transcript_id", "_", type, "_", pred$start, ":", pred$end )
  mcols(tr_copy)$"transcript_name" <- paste0( mcols(tr_copy)$"transcript_name", "_", type, "_", pred$start, ":", pred$end )

  copy <- tr_copy[  match(paste0( tr_ids, "_", type, "_", pred$start, ":", pred$end ),  mcols(tr_copy)$transcript_id )]
  ranges(copy) <- IRanges(pred$start, pred$end)  ## replace the ranges
  mcols(copy)[c("exon_number",  "exon_version")] <- NA  ## remove the exon id, number and version
  mcols(copy)$"exon_id" <- paste0( mcols(copy)$"exon_id", "_", type, "_", pred$start, ":", pred$end )

  return(c(tr_copy, copy) )
}

novel_entries <- lapply(seq_along(tr_ids), function(x) if(length(tr_ids[[x]])>0){ new_transcript(exon_pred[x],tr_ids[[x]], exons = exons) } else NA )


## remove the novel exons that do not have a transcript
novel_entries <- do.call("c", novel_entries[!is.na(novel_entries)])

### add the novel entries to the original gtf file
gtf_added <- c(gtf, novel_entries)

export(gtf_added, OUTFILE)

