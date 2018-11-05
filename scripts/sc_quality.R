## Compare the Phred quality of the first soft-clipped base in each read in a BAM file.


## TODO: what about strandedness? is the quality in the same orientation than the 
## CIGAR?
## or does this depend on the strand?


#' Quality scores at index
#' 
#' Extract the quality scores at the specified index for soft-clipped regions at
#' the 5' of 3' end of a read
#'
#' @param qualities character vector with Phred qualities in ASCII format
#' @param location string either `five_prime` or `three_prime`
#' @param i integer desired index in the quality strings 
#'
#' @return integer vector with the quality scores at index i in all quality strings
#' @export
#'
#' @examples
qualities_at <- function(qualities, i, location){
  try(if (!location %in% c("five_prime", "three_prime")) 
    stop("read must be either `five_prime` or `three_prime`"))
  
  if(location == "five_prime"){
    pos <- ifelse(nchar(qualities)-i>0, 
                  nchar(qualities)-i, 
                  0)
  } else {
    pos <- ifelse(i+1<=nchar(qualities), 
                  i+1, 
                  0)
  }
  q <- substr(qualities, pos, pos)
  unlist(as(PhredQuality(q), "IntegerList"))
}


#' Function to get the maximum and the number of elements from an input vector
#'
#' @param d numeric vector
#'
#' @return named vector with "y" = maximum of d and "label" = number of elements
#'  in d
#' @export
#'
#' @examples
give_n <- function(d){
  return(c(y = max(d)*1.05, label = length(d))) 
  # experiment with the multiplier to find the perfect position
}

#' Plot quality scores
#' 
#' Boxplot and violin plot of the quality scores at different positions relative 
#' to the start of soft-clipped regions.
#'
#' @param bam path to BAM file
#' @param outprefix string filepath and prefix of the output plot
#'
#' @return 
#' @export
#'
#' @examples
plot_quality_scores <- function(bam, outprefix, title = ""){
  aln <- readGAlignmentPairs(bam, strandMode=2, use.names=TRUE, 
                             param = ScanBamParam(what = "qual" ))
  
  ## Only keep all soft-clipped reads
  s1 <- GenomicAlignments::first(aln)[grepl("S", cigar(GenomicAlignments::first(aln)))]
  s2 <- GenomicAlignments::second(aln)[grepl("S", cigar(GenomicAlignments::second(aln)))]
  
  ## Find the location of the "S" in the read: https://support.bioconductor.org/p/75307/
  s1_r <- cigarRangesAlongQuerySpace(cigar(s1), ops="S")
  s2_r <- cigarRangesAlongQuerySpace(cigar(s2), ops="S")
  
  ## We split the soft-clipped ranges into 5' and 3'
  s1_r_5 <- s1_r[start(s1_r) == 1]
  s1_r_3 <- s1_r[start(s1_r) != 1]
  s2_r_5 <- s2_r[start(s2_r) == 1]
  s2_r_3 <- s2_r[start(s2_r) != 1]
  
  ## Histogram of the length of the soft-clipped regions =======================
  s1_width <- c(unlist(width(s1_r_3)), unlist(width(s1_r_5)))
  s2_width <- c(unlist(width(s2_r_3)), unlist(width(s2_r_5)))
  
  dat_width <- data.frame(sc_length = c(s1_width, s2_width), 
                          read = c(rep("first", length(s1_width)), 
                                   rep("second", length(s2_width))))
  
  p <- ggplot(dat_width, aes(x=sc_length, fill=read)) +
    geom_histogram(position="dodge", binwidth = 1) +
    theme_bw() +
    xlab("length of the soft-clipped region") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(outprefix, "_sc_length_histogram.pdf"), p, 
         width=7, height = 5)
  
  ## ===========================================================================

  ## expand by one position to the right or left so we can include the quality 
  ## values of the base preceding the soft-clipped region
  s1_r_5 <- resize(s1_r_5, width = width(s1_r_5)+1, fix = "start")
  s1_r_3 <- resize(s1_r_3, width = width(s1_r_3)+1, fix = "end")
  s2_r_5 <- resize(s2_r_5, width = width(s2_r_5)+1, fix = "start")
  s2_r_3 <- resize(s2_r_3, width = width(s2_r_3)+1, fix = "end")
  s1_qual_5 <- as.character(mcols(s1)$qual[s1_r_5])
  s1_qual_3 <- as.character(mcols(s1)$qual[s1_r_3])
  s2_qual_5 <- as.character(mcols(s2)$qual[s2_r_5])
  s2_qual_3 <- as.character(mcols(s2)$qual[s2_r_3])
  
  # s1_qual_5 <- as(mcols(s1)$qual[s1_r_5], "IntegerList")
  # s1_qual_3 <-as(mcols(s1)$qual[s1_r_3], "IntegerList")
  # s1_qual_5 <- s1_qual_5[sapply(s1_qual_5, function(x) length(x) !=0)]
  # s1_qual_3 <- s1_qual_3[sapply(s1_qual_3, function(x) length(x) !=0)]
  s1_qual_5 <- s1_qual_5[nchar(s1_qual_5)>0]
  s1_qual_3 <- s1_qual_3 [nchar(s1_qual_3)>0]
  s2_qual_5 <- s2_qual_5[nchar(s2_qual_5)>0]
  s2_qual_3 <- s2_qual_3 [nchar(s2_qual_3)>0]
  
  ## Quality values per soft-clipped region at different idices
  ## the base at index 0 is NOT soft-clipped 
  dat <- list()
  for (i in c(0, 1, 2, 3)){
    q <- qualities_at(s1_qual_5, i, "five_prime")
    dat[["qual"]] <- c(dat[["qual"]], q)
    dat[["position"]] <- c(dat[["position"]], rep(i, length(q)))
    dat[["read"]] <- c(dat[["read"]], rep("first", length(q)))
    
    q <- qualities_at(s1_qual_3, i, "three_prime")
    dat[["qual"]] <- c(dat[["qual"]], q)
    dat[["position"]] <- c(dat[["position"]], rep(i, length(q)))
    dat[["read"]] <- c(dat[["read"]], rep("first", length(q)))
    
    q <- qualities_at(s2_qual_5, i, "five_prime")
    dat[["qual"]] <- c(dat[["qual"]], q)
    dat[["position"]] <- c(dat[["position"]], rep(i, length(q)))
    dat[["read"]] <- c(dat[["read"]], rep("second", length(q)))
    
    q <- qualities_at(s2_qual_3, i, "three_prime")
    dat[["qual"]] <- c(dat[["qual"]], q)
    dat[["position"]] <- c(dat[["position"]], rep(i, length(q)))
    dat[["read"]] <- c(dat[["read"]], rep("second", length(q)))
  }
  dat <- data.frame(quality_score = dat[["qual"]], 
                    position = as.factor(dat[["position"]]), 
                    read = as.factor(dat[["read"]]))
  
  dodge <- position_dodge(width = 0.8)
  p <- ggplot(dat, aes(x = position, y = quality_score, color = read)) +
    geom_boxplot(position=dodge, width=0.3, outlier.shape = NA) +
    geom_violin(position=dodge, scale="area", alpha =0.1) +
    theme_bw() +
    xlab("position relative to first soft-clipped base") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(fun.data = give_n, geom = "text", position = dodge)  
  ggsave(paste0(outprefix, "_quality_scores_per_position.pdf"), p, 
         width=7, height = 6)

  ## TODO: write table with percentage of soft-clipped bases in the sample
}


library(GenomicAlignments)
library(dplyr)
library(ggplot2)

## parameters
BAM <- snakemake@input[["bam"]]
OUTPREFIX <- snakemake@params[["outprefix"]]
TITLE <- snakemake@params[["title"]]

# BAM <- "../simulation/mapping/STAR/me_exon/default/pass2_Aligned.out_s.bam"
# OUTPREFIX <- "../simulation/mapped_truth/star/me_exon/default/star_"
# TITLE <- "STAR: me_exon, default"

plot_quality_scores(BAM, OUTPREFIX, TITLE)
