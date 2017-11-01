#!/usr/local/R/R-3.4.0/bin/Rscript

#SERVER <- '/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch,user=vbarbo'      # if running in my computer
# SERVER <- NULL                                                                 # if running in the server
# GTF_FILE <- '/home/Shared_sherborne/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/include_mic_full_coordinate_gtf/Homo_sapiens.GRCh37.85_chr19_22_me_included_all_trans.gtf'
# SIM_ISOFORMS_RESULTS <- '/home/Shared_sherborne/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/simulation/SRR3192428/simulated_data/chr19_22_SRR3192428.sim.isoforms.results'
# FASTQ_FILE_READ1 <- '/home/Shared_sherborne/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/simulation/SRR3192428/simulated_data/chr19_22_SRR3192428_1.fq'
# FASTQ_FILE_READ2 <- '/home/Shared_sherborne/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/simulation/SRR3192428/simulated_data/chr19_22_SRR3192428_2.fq'
# OVERLAP <- 1
# MAX_MIC_LENGTH <- Inf
# OUTPUT_FILE <- '/home/Shared_sherborne/data/seq/microexon_simulation/microexon_study/new_simulation/script_output/exons_truth/exons_truth_all_exons.txt'
# READS_LENGTH <- 101
# CORES<-25



print(version)

GTF_FILE <- snakemake@input[["gtf"]]
SIM_ISOFORMS_RESULTS <- snakemake@input[["sim_iso_res"]]
FASTQ_FILE_READ1 <- snakemake@input[["fastq1"]]
FASTQ_FILE_READ2 <- snakemake@input[["fastq2"]]
OUTPUT_FILE <- snakemake@output[[1]]
CORES<- snakemake@threads

OVERLAP <- 1
MAX_MIC_LENGTH <- Inf
READS_LENGTH <- 101


library(parallel)

### sim.iso.res = sim.isoforms.results file.
sim.iso.res <- read.table(SIM_ISOFORMS_RESULTS,h=T)
sim.iso.res <- cbind(sim.iso.res,sid=1:nrow(sim.iso.res))


### GTF file with sid
library(rtracklayer)
gtf <- import(GTF_FILE)
gtf <- gtf[mcols(gtf)$type=='exon']
mcols(gtf)$sid <- merge(as.data.frame(gtf),sim.iso.res[,c('transcript_id','sid')],by='transcript_id',sort=F)$sid   # Adding the matecolumn sid


### The sequence of the read 1.
# dir=0, gtf=+ --> normal
# dir=0, gtf=- --> reversed complementrary
# dir=1, gtf=+ --> reversed complementrary
# dir=1, gtf=- --> normal



inds <- split(1:length(gtf),gtf$sid)
inds1 <- sapply(inds,.subset,1)
gtfs <- split(gtf, mcols(gtf)$sid)
ws <- width(gtfs)
str <- as.character(strand(gtf))[inds1]
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
mcols(gtf)$pos_start <- unsplit(lapply(ps,function(u) u[,1]), gtf$sid)
mcols(gtf)$pos_end <- unsplit(lapply(ps,function(u) u[,2]), gtf$sid)

rm(gtfs,inds,inds1,ps,s,str,getPositions)


### output = FASTQ file with positions for read1 and read2 
library("ShortRead")
fq <- readFastq(FASTQ_FILE_READ1)
output1 <- as.vector(id(fq))
output1 <- sapply(output1, function(x) unlist(strsplit(x,'/'))[1])
output1 <- t(sapply(output1, function(x) unlist(strsplit(x,'_'))))
mode(output1) <- 'numeric'
output1 <- cbind(output1,as.data.frame(sread(fq)))

fq <- readFastq(FASTQ_FILE_READ2)
output2 <- as.vector(id(fq))
output2 <- sapply(output2, function(x) unlist(strsplit(x,'/'))[1])
output2 <- t(sapply(output2, function(x) unlist(strsplit(x,'_'))))
mode(output2) <- 'numeric'
output2 <- cbind(output2,as.data.frame(sread(fq)))
colnames(output2) <- colnames(output1) <- c('rid','dir','sid','pos','insertL','seq')
rownames(output2) <- rownames(output1) <- NULL
output1 <- output1[output1$sid!=0,]
output2 <- output2[output2$sid!=0,]
rm(fq)

output <- output1[,1:5]
names(output)[5] <- 'insert'
strand <- split(as.vector(strand(gtf)), gtf$sid)
strand <- strand[names(strand) %in% output$sid]
strand <- lapply(strand, unique)
strand <- unsplit(strand, output$sid)
output <- cbind(output,gtf=strand)



pos <- data.frame(start1=numeric(nrow(output)), end1=numeric(nrow(output)), start2=numeric(nrow(output)), end2=numeric(nrow(output)))
foward <- (output$dir==0 & output$gtf=='+') | (output$dir==1 & output$gtf=='-')

pos[foward, 'start1'] <- output$pos[foward] + 1
pos[foward, 'end1'] <- pos[foward, 'start1'] + READS_LENGTH - 1
pos[foward, 'end2'] <- pos[foward, 'start1'] + output$insert[foward] - 1
pos[foward, 'start2'] <- pos[foward, 'end2'] - READS_LENGTH + 1

trans_length <- lapply(ws, sum)
trans_length <- trans_length[names(trans_length) %in% output$sid]
trans_length <- unsplit(trans_length,output$sid)

pos[!foward, 'end1'] <- trans_length[!foward] - output$pos[!foward]
pos[!foward, 'start1'] <- pos[!foward, 'end1'] - READS_LENGTH + 1
pos[!foward, 'start2'] <- pos[!foward, 'end1'] - output$insert[!foward] + 1
pos[!foward, 'end2'] <- pos[!foward, 'start2'] + READS_LENGTH - 1

output <- cbind(output,pos)


### Count of reads for each microexon (lines in the GTF file)
gtf_chr <- gtf[width(gtf)<=MAX_MIC_LENGTH]

sids <- split(1:length(output$sid), output$sid)
cr <- mclapply(mc.cores=CORES, 1:length(gtf_chr), function(i){
  s <- which(names(sids) == gtf_chr$sid[i])
  if(length(s)>0){
    cr <- table((gtf_chr$pos_end[i] - output$start1[sids[[s]] ] >= (OVERLAP-1) & output$end1[sids[[s]] ] - gtf_chr$pos_start[i] >= (OVERLAP-1)) |
                  (gtf_chr$pos_end[i]-output$start2[sids[[s]] ] >= (OVERLAP-1) & output$end2[sids[[s]] ] - gtf_chr$pos_start[i] >= (OVERLAP-1)))['TRUE']   # All reads that contain the microexons i or at least OVERLAP nucleotides of it.
    if(is.na(cr)) cr <- 0
    cr
  }else{
    0
  }
})
cr<-sapply(cr,function(x) setNames(x,NULL))
mcols(gtf_chr)$count_reads<-cr

exons_count <- subset(as.data.frame(gtf_chr),select=c(seqnames,start,end,exon_id,transcript_id,sid,count_reads))


### Count of reads for each microexon (locations)
loc <- unique(exons_count[,1:3])
loc$length <- loc$end - loc$start +1
count_reads <- mclapply(mc.cores=CORES, 1:nrow(loc), function(i) sum(exons_count$count_reads[exons_count$seqnames==loc$seqnames[i] & exons_count$start==loc$start[i] & exons_count$end==loc$end[i]]))
count_reads <- sapply(count_reads, function(x) x)
loc <- cbind(loc,count_reads)


### Microexons that are present in the data
# loc_present <- loc[loc$count_reads>0,]
# loc_present <- loc_present[order(loc_present$seqnames,loc_present$start,loc_present$end),]
# loc_present <- cbind(nr=1:nrow(loc_present),loc_present)


### Writing
#write.table(loc_present,paste0(SERVER,OUTPUT_FILE),quote=F,row.names=F,sep='\t')


# also microexons with expression = 0
loc <- loc[order(loc$seqnames,loc$start,loc$end),]
loc <- cbind(nr=1:nrow(loc),loc)
write.table(loc,OUTPUT_FILE,quote=FALSE,row.names=FALSE,sep='\t')
