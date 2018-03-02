### This script takes a GTF and genome fasta file as input and writes the transcript sequences in fasta format.

GTF <- snakemake@input[["gtf"]]
GENOME <- snakemake@input[["genome"]]
OUTFILE <- snakemake@output[["outfile"]]

# GENOME <- "/home/Shared/data/seq/microexon_simulation/GRCh37.82/genome/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
# GTF <- "/home/Shared/kathi/microexon_pipeline/simulation/reduced_GTF_with_predicted_exons/me_exon/GRCh37.85_chr19_22_novel_exons_outSJfilterOverhangMin6.gtf"
# OUTFILE <- "/home/Shared/kathi/microexon_pipeline/simulation/transcriptome/test/test.fasta"
# "me_exon/GRCh37.85_chr19_22_novel_exons_outSJfilterOverhangMin6.fasta"


## my solution, faster if we only have transcripts on a few chromosomes (e.g. the simulation)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)

# seq_kathi <- function(){
fa <- open(FaFile(GENOME))
idx <- scanFaIndex(fa) ## read the fasta index for fast accession of single chromosomes
gtf <- import(GTF) 

exons <- gtf[ mcols(gtf)$type == "exon" ]
exons_seqnames <- split(exons, seqnames(exons)) ## split all exons according to seqnames

### For a given chromosome and strand, extract all transcript sequences
## chr_seq: DNAStringSet with chromosome sequence
## stra: strand, either "+" or "-"
## exons: GRanges object with all exons on the chromosome
get_trans_seq_stranded <- function(chr_seq, stra, exons){
	ex <- exons_seqnames[[chr]]
	ex <- ex[ strand(ex) == stra ] 
	ex <-  data.frame(start=start(ex),end=end(ex),transcript_id=mcols(ex)$transcript_id)
	if(stra == "+"){
		ex <-  ex[order(ex$start),] ## order ascending
	} else{  ## "-" strand
		ex <-  ex[order(ex$start, decreasing=TRUE),] ## order decreasing
	}
	tr <- split(IRanges(start = ex$start, end = ex$end), ex$transcript_id)
	extractTranscriptSeqs(chr_seq, tr, strand=stra)
}

tr_seqs <- DNAStringSet()
## go through all chromosomes in the gtf file and extract the transcript sequences
for(chr in names(exons_seqnames)){
	ex <- exons_seqnames[[chr]]
	chr_seq <- unlist(getSeq(fa, param=idx[ which(seqnames(idx) == chr) ]) )
	tr_seqs <- c(tr_seqs, get_trans_seq_stranded(chr_seq, "+", ex))
	tr_seqs <- c(tr_seqs, get_trans_seq_stranded(chr_seq, "-", ex))
}
# }

writeXStringSet(tr_seqs, OUTFILE, format="fasta")







### Marks solutions, might be faster if we have transcripts on most chromosomes
# seq_mark <- function() {
# 	ref <- readDNAStringSet(GENOME)
# 	names(ref) <- unlist(lapply( strsplit(names(ref), " "), function(x) x[1]))

# 	## order the exons on the + and - strand
# 	gtf <- import(GTF) 
# 	exons <- gtf[ mcols(gtf)$type == "exon" ]
# 	exons_plus <- exons[strand(exons) == "+"]
# 	exons_plus <- sort(exons_plus)
# 	exons_neg <- exons[strand(exons) == "-"]
# 	exons_neg <- sort(exons_neg, decreasing=TRUE)
# 	exons <- c(exons_plus, exons_neg)

# 	regs <- split(exons, mcols(exons)$transcript_id)

# 	seqs <- getSeq(ref, regs)
# 	tr_seqs <- DNAStringSet(sapply(seqs, paste, collapse=""))
# }

## we read in the whole genome fasta in Marks method, but in my method we only read the chromosomes that actually have genes in the gtf file.
## This is one of the reasons why my method is much faster


## much faster:
### using DNAStringSet and GRanges as input to extractTranscriptSeq
# seq_mark2 <- function(){
# 	library(rtracklayer) ## for import()
# 	library(Biostrings) ## for readDNAStringSet()
# 	library(GenomicFeatures) ## for extractTranscriptSeqs()
# 	library(BSgenome) ## for DNAStringSet input to extractTranscriptSeqs()
# 	gtf <- import(GTF) 
# 	exons <- gtf[ mcols(gtf)$type == "exon" ]
# 	exons_plus <- exons[strand(exons) == "+"]
# 	exons_plus <- sort(exons_plus)
# 	exons_neg <- exons[strand(exons) == "-"]
# 	exons_neg <- sort(exons_neg, decreasing=TRUE)
# 	exons <- c(exons_plus, exons_neg)

# 	annos <- split(exons, mcols(exons)$transcript_id)

# 	ref <- readDNAStringSet(GENOME)
# 	names(ref) <- unlist(lapply( strsplit(names(ref), " "), function(x) x[1]))
# 	tr_seqs <- extractTranscriptSeqs(ref, annos)
# }

# > microbenchmark( seq_mark(), seq_mark2(), seq_kathi(), times = 3 )
# Unit: seconds
#         expr       min        lq     mean    median        uq       max neval
#   seq_mark() 67.719713 68.833258 70.51227 69.946803 71.908551 73.870300     3
#  seq_mark2() 33.583519 34.732391 35.27168 35.881263 36.115756 36.350248     3
#  seq_kathi()  8.028004  8.120758  8.43014  8.213513  8.631208  9.048903     3



# ## solution with dplyr
# mdata = tibble(start=start(ex_neg),end=end(ex_neg),tid=mcols(ex_neg)$transcript_id)
# msorted = mdata %>% arrange(start,.by_group=T)