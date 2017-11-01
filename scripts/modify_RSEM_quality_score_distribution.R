## This script reads the quality score distribution for an RSEM .model file and changes the pobability for a quality score of 2 to 0 to prevent sampling reads with overall low quality
## 20.05.2016

markov_prob <- snakemake@input[[1]]
OUTFILE <- snakemake@output[[1]]

d <- read.table(markov_prob)
d[,3] <- rep(0,length(d[,3]))  # set the initial probability for quality score 2 to 0: and set the transition to and from 2 to 0
# renormalize the probabilities per row to 1
sumvector <- apply(d, 1, sum)  # get the sum of all rows
d <- d/sumvector  # divide the probability matrix by the sums per row to get the renormalized probabilities
d[is.na(d)] <- 0

write.table(d, file = OUTFILE, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
