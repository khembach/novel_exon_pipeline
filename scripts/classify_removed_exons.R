## This script classifies the removed exons into 
## "easy": novel cassette exon with 2 novel splice junctions
## "medium": novel cassette exons with only 1 novel splice junction
## "complicated": novel exon without a novel splice junction, because its ends overlap with annotated exons

#### TOD what about the terminal exons? if their only junction is already annotated they should be classified as complicated and not medium!!!
## or should all of them be medium??


ME_EXON <- "simulation/reduced_GTF/removed_microexons_exons_unique.txt"
# ME <- ""
# EXON <- ""

removed <- read.table(ME_EXON, header=TRUE) ## 261 removed exons

removed$type <- ifelse(removed$width <= 27, "me", "exon")
removed$class <- ifelse(removed$unique_exon, "easy", NA)
removed$class[which( removed$shared_exon_start >0 & removed$shared_exon_end >0 )] <- "complicated"
removed$class[(removed$shared_exon_start >0 & removed$shared_exon_end == 0  | removed$shared_exon_start == 0 & removed$shared_exon_end > 0  )] <- "medium"
removed$class[ which( is.na( removed$class ) ) ] <- "medium"   ## exons that have unique start and end, but still overlap with another exon (most of them are rather easy to find)

### terminal exons:
terminal <- which(is.na(removed$lend) | is.na(removed$rstart) )
removed$terminal <- FALSE
removed$terminal[terminal] <- TRUE

removed$class[terminal] <- "medium"
removed$class[ terminal[  which(removed$shared_exon_start[terminal] >0 | removed$shared_exon_end[terminal] >0 ) ] ] <- "complicated"





# > table(removed$class)
# complicated        easy      medium
#          59          94         108


me <- removed[removed$width <= 27,]
exon <- removed[removed$width > 27,]
# > table(me$class)
##me
# complicated        easy      medium
#          47          38          39

## exon
# complicated        easy      medium
#          12          56          69

## only take the exon location into account, not to exons it connects to in the transcript (should be 100me and 100 exons)

removed_unique <- removed[,-which(colnames(removed) %in% c("rstart", "lend"))]

## remove all duplicate rows 
removed_unique <- removed_unique[ rownames(unique(removed_unique[, c("seqnames", "start", "end", "strand")])),]

me_unique <- removed_unique[removed_unique$width <= 27,]
exon_unique <- removed_unique[removed_unique$width > 27,]

# table(me_unique$class)
# table(exon_unique$class)
# ## me
# complicated        easy      medium
#           3          42          49
# ## exon
# complicated        easy      medium
#           8          49          43


write.table(removed, "simulation/reduced_GTF/removed_microexons_exons_unique_classified.txt", sep = "\t", row.names = FALSE, quote=FALSE )