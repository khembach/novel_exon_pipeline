### Plot precision recall curves for the novel exons.
## First, we evaluate the exons by class. Then we only evaluate the exons with at least 1 mapped junction read (all others are not expected to be predicted). 
## Last, we we separate the exons according to their location, type or expression.


library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)


REMOVED <- snakemake@input[["removed"]]  ## al removed exons, with mapped junction counts
PREDICTION <- snakemake@input[["prediction"]]
OUTDIR <- snakemake@output[["outdir"]]


## stringtie
# REMOVED = "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/mapped_junction_count/removed_me_exon_unique_classified_outSJfilterOverhangMin6_junc_count.txt"
# PREDICTION = "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/stringtie/me_exon/novel_exons_outSJfilterOverhangMin6_stringtie.txt"



# REMOVED <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/mapped_junction_count/removed_me_exon_unique_classified_outSJfilterOverhangMin6_junc_count.txt"
# REMOVED <- "simulation/reduced_GTF/removed_microexons_exons_unique_classified.txt"
# PREDICTION <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/filtered_SJ/me_exon/novel_exons_outSJfilterOverhangMin6.txt"
# PREDICTION <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/filtered_SJ/two_junc_reads/me_exon/novel_exons_outSJfilterOverhangMin6.txt"

# params <- c("default", "outSJfilterOverhangMin9", "outSJfilterOverhangMin6", "outSJfilterCountTotalMin3", "scoreGenomicLengthLog2scale0", "alignSJoverhangMin3")
# PRED_PATH <- "simulation/analysis/filtered_SJ/me_exon/novel_exons_"
# OUTDIR <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/analysis/exon_prediction_performance/PR/test/"


# SJFILE <- "/Volumes/Shared/kathi/microexon_pipeline/simulation/mapping/STAR/me_exon/outSJfilterOverhangMin6/pass2_SJ.out.tab"
# sj <- fread(SJFILE)
# colnames(sj) <- c("seqnames", "start", "end", "strand", "motif", "annotated", "unique", "mutimapping", "maxoverhang")
# 
# #strand: (0: undefined, 1: +, 2: -)
# sj$strand[sj$strand==0] <- "*"
# sj$strand[sj$strand==1] <- "+"
# sj$strand[sj$strand==2] <- "-"




novelExons <- read.table(PREDICTION, header=TRUE, stringsAsFactors = FALSE)
##  list of all removed exons (the truth)
removed <- read.table(REMOVED, header=TRUE) ## 261 removed exons

######## plot Precision-Recall curve stratified by the min # of junction reads
## label all predictions/truth: if cassette -> 4 coordinates, if right terminal: 2 coordinates; if left terminal: 2 coordinates
novelExons$location <- "cassette"
novelExons$location[is.na( novelExons$lend )] <- "left_terminal"
novelExons$location[is.na (novelExons$rstart)] <- "right_terminal"

removed$location <- "cassette"
removed$location[is.na( removed$lend )] <- "left_terminal"
removed$location[is.na (removed$rstart)] <- "right_terminal"
removed$ID <- 1:nrow(removed)

evaluate_prediction <-function(criteria, pred, truth){
  ## check if the predictions can be found among the truth and add a column "label" with TRUE, otherwise NA
  labeled <- left_join( pred %>% select(criteria),  truth %>% select(c(criteria, "ID")) %>% mutate(label = TRUE) ) %>% select(matched_ID = ID, label)
  labeled$label[ is.na(labeled$label) ] <- FALSE
  cbind(pred, labeled)
}


###### add class and true junc counts, true gene count and add the matched truth ID
## label all predictions with TRUE or FALSE depending if the prediction is correct or not
cassette <- c("seqnames", "start", "end", "strand", "lend", "rstart")
left_terminal <- c("seqnames",  "end", "strand", "lend", "rstart")
right_terminal <- c("seqnames",  "start", "strand", "lend", "rstart")
eval_criteria <- list(cassette = cassette, left_terminal = left_terminal, right_terminal = right_terminal)

novelExons_split <- split(novelExons, novelExons$location)
novelExons_split <- lapply(seq_along(eval_criteria), function(x) evaluate_prediction(eval_criteria[[x]], novelExons_split[[names(eval_criteria)[x]]], removed))
novelExons <- unsplit(novelExons_split,  novelExons$location)
  
  
## given a table with labeled predicted exons, compute precision and recall
compute_pr_rec <- function(novelExons, nr_removed){
  ## compute TP and FP for each min_read bin
  eval_dat <- novelExons %>%group_by(min_reads) %>%
    summarise(tp = sum(label), fp = sum(!label)) %>% arrange(desc(min_reads)) %>%
    mutate(sum_tp = cumsum(tp), sum_fp = cumsum(fp)) 
  ## compute FN
  eval_dat <- eval_dat %>% mutate(fn = nr_removed - sum_tp) 
  ## compute PR and REC
  eval_dat <- eval_dat %>% mutate(precision = sum_tp/(sum_tp+sum_fp), recall = sum_tp/(sum_tp + fn))
  eval_dat
}

####### all predictions, all removed exons
dat_all <- compute_pr_rec(novelExons, nrow(removed))
precision_all <- dat_all %>% select(min_reads, precision)

p <- ggplot(dat_all, aes(x=recall, y = precision)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + xlim(c(0,1)) + ggtitle(paste0("All removed exons (", nrow(removed), ")"))
ggsave(paste0(OUTDIR, "PR.png"), p, device = "png") 

## PR vs reads
p <- ggplot(dat_all, aes(x = min_reads, y = precision)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + ggtitle(paste0("All removed exons (", nrow(removed), ")"))
ggsave(paste0(OUTDIR, "pr_reads.png"), p, device = "png")

## PR vs reads
p <- ggplot(dat_all, aes(x = min_reads, y = recall)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + ggtitle(paste0("All removed exons (", nrow(removed), ")"))
ggsave(paste0(OUTDIR, "rec_reads.png"), p, device = "png")


######### class: easy, medium, hard
id_class <- split(removed$ID, removed$class)  ## all IDs that belong to a certain class

dat_class <- list()
for(cl in names(id_class)){
  dat_class[[cl]] <- novelExons %>% filter(matched_ID %in% id_class[[cl]])
}

## given a table with labeled predicted exons, compute precision and recall
compute_rec <- function(novelExons, nr_removed){
  ## compute TP and FP for each min_read bin
  eval_dat <- novelExons %>%group_by(min_reads) %>%
    summarise(tp = sum(label)) %>% arrange(desc(min_reads)) %>%
    mutate(sum_tp = cumsum(tp)) 
  ## compute FN
  eval_dat <- eval_dat %>% mutate(fn = nr_removed - sum_tp) 
  ## REC
  eval_dat <- eval_dat %>% mutate(recall = sum_tp/(sum_tp + fn))
  eval_dat
}

dat_class <- lapply( names(id_class), function(x) compute_rec( dat_class[[x]], length(id_class[[x]]) ) ) ## compute recall for all classes
names(dat_class) <- names(id_class)
dat_class <- rbindlist(dat_class, idcol = "class")
dat_class <- left_join(dat_class, precision_all) ## add the corresponding precision value from all predictions 

p <-  ggplot(dat_class, aes(x = min_reads, y =recall, color = class)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + geom_point() + ggtitle(paste0("All removed exons (", nrow(removed), ")")) 
ggsave(paste0(OUTDIR, "rec_reads_class.png"), p, device = "png")

p <- ggplot(dat_class, aes(x = recall, y = precision, color = class))  + geom_path(size = 1.5) + theme_bw() + geom_point()+ ylim(c(0, 1)) + xlim(c(0,1)) + ggtitle(paste0("All removed exons (", nrow(removed), ")"))
ggsave(paste0(OUTDIR, "PR_class.png"), p, device = "png")



######## We expact that the missed exons are all exons that do not have mapped junction reads and are thus impossible to predict.
## separate all exons according to mapped number of reads crossing the SJ: 0 reads or >0 junction reads

## split the exons according to the minimal number of mapped reads over each junction (0 reads or more)
id_min_junc <- split(removed$ID, pmin(removed$count_lsj, removed$count_rsj, na.rm = TRUE) == 0)
names(id_min_junc) <- ifelse(names(id_min_junc) == "TRUE", "no_junction_reads", "junction_reads")

dat_min_junc <- list()
for(cl in names(id_min_junc)){
  dat_min_junc[[cl]] <- novelExons %>% filter(matched_ID %in% id_min_junc[[cl]])
}

dat_min_junc <- lapply( names(id_min_junc), function(x) compute_rec( dat_min_junc[[x]], length(id_min_junc[[x]]) ) ) ## compute recall for all classes
names(dat_min_junc) <- names(id_min_junc)
dat_min_junc <- rbindlist(dat_min_junc, idcol = "after_mapping")
dat_min_junc <- left_join(dat_min_junc, precision_all) ## add the corresponding precision value from all predictions 

p <-  ggplot(dat_min_junc, aes(x = min_reads, y =recall, color = after_mapping)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1))+ ggtitle("Exons with >0 mapped junction reads")
ggsave(paste0(OUTDIR, "rec_junc_read.png"), p, device = "png")

p <- ggplot(dat_min_junc, aes(x = recall, y = precision, color = after_mapping))  + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + xlim(c(0,1))  + geom_point() + ggtitle("All removed exons")+ ylab("precision (all exons)")
ggsave(paste0(OUTDIR, "PR_junc_reads.png"), p, device = "png")


####### From here on, we only consider the exons with at least 1 junction read
## compute the precision for only the exons with junction reads
removed <- removed %>% filter(pmin(removed$count_lsj, removed$count_rsj, na.rm = TRUE) > 0)  ## all exons with at least 1 junction read
dat_all <- compute_pr_rec(novelExons, nrow(removed))
precision_all <- dat_all %>% select(min_reads, precision)


id_class <- split(removed$ID, removed$class)  ## all IDs that belong to a certain class

dat_class <- list()
for(cl in names(id_class)){
  dat_class[[cl]] <- novelExons %>% filter(matched_ID %in% id_class[[cl]])
}

dat_class <- lapply( names(id_class), function(x) compute_rec( dat_class[[x]], length(id_class[[x]]) ) ) ## compute recall for all classes
names(dat_class) <- names(id_class)
dat_class <- rbindlist(dat_class, idcol = "class")
dat_class <- left_join(dat_class, precision_all) ## add the corresponding precision value from all predictions 

p <-  ggplot(dat_class, aes(x = min_reads, y =recall, color = class)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + geom_point()+ ggtitle("Exons with >0 mapped junction reads")
ggsave(paste0(OUTDIR, "rec_reads_class_expr.png"), p, device = "png")

p <- ggplot(dat_class, aes(x = recall, y = precision, color = class))  + geom_path(size = 1.5) + theme_bw() + geom_point()+ ylim(c(0, 1)) + xlim(c(0,1))+ ggtitle("Exons with >0 mapped junction reads")+ ylab("precision (all mapped exons)")
ggsave(paste0(OUTDIR, "PR_class_expr.png"), p, device = "png")



###### me, normal exons
id_type <- split(removed$ID, removed$type)  ## all IDs that belong to a certain class
novelExons <- novelExons %>% mutate(type=ifelse((end-start+1)<=27, "me", "exon"))

dat_type <- list()
for(ty in names(id_type)){
  dat_type[[ty]] <- novelExons %>% filter(type == ty)
}

dat_type <- lapply( names(id_type), function(x) compute_pr_rec( dat_type[[x]], length(id_type[[x]]) ) ) ## compute recall for all classes
names(dat_type) <- names(id_type)
dat_type <- rbindlist(dat_type, idcol = "type")

p <-  ggplot(dat_type, aes(x = min_reads, y =recall, color = type)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1))+ ggtitle("Exons with >0 mapped junction reads")
ggsave(paste0(OUTDIR, "rec_reads_exon_type.png"), p, device = "png")

p <- ggplot(dat_type, aes(x = recall, y = precision, color = type ))  + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + xlim(c(0,1)) + ggtitle("Exons with >0 mapped junction reads")+ ylab("precision (all mapped exons)")
ggsave(paste0(OUTDIR, "PR_exon_type.png"), p, device = "png")



###### cassette exons, terminal exon
# id_location <- split(removed$ID, removed$location)  ## all IDs that belong to a certain class
id_location <- list()
id_location[["cassette"]] <- removed$ID[removed$location == "cassette"]
id_location[["terminal"]] <- removed$ID[removed$location != "cassette"]
  
dat_location <- list()
dat_location[["cassette"]] <- novelExons %>% filter(location == "cassette")
dat_location[["terminal"]] <- novelExons %>% filter(location != "cassette")

dat_location <- lapply( names(id_location), function(x) compute_pr_rec( dat_location[[x]], length(id_location[[x]]) ) ) ## compute recall for all classes
names(dat_location) <- names(id_location)
dat_location <- rbindlist(dat_location, idcol = "location")

p <-  ggplot(dat_location, aes(x = min_reads, y =recall, color = location)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + ggtitle("Exons with >0 mapped junction reads")
ggsave(paste0(OUTDIR, "rec_reads_exon_location.png"), p, device = "png")

p <- ggplot(dat_location, aes(x = recall, y = precision, color = location ))  + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + xlim(c(0,1)) + ggtitle("Exons with >0 mapped junction reads")+ ylab("precision (all mapped exons)")
ggsave(paste0(OUTDIR, "PR_exon_location.png"), p, device = "png")



###### low and high expressed exons: simulated number of reads below or above median
median_expression <- median(removed$count_reads)
id_expr <- split(removed$ID, removed$count_reads< median_expression )  ## all IDs that belong to a certain class
names(id_expr) <- c("high", "low")

dat_expr <- list()
for(cl in names(id_expr)){
  dat_expr[[cl]] <- novelExons %>% filter(matched_ID %in% id_expr[[cl]])
}

dat_expr <- lapply( names(id_expr), function(x) compute_rec( dat_expr[[x]], length(id_expr[[x]]) ) ) ## compute recall for all classes
names(dat_expr) <- names(id_expr)
dat_expr <- rbindlist(dat_expr, idcol = "expression")
dat_expr <- left_join(dat_expr, precision_all) ## add the corresponding precision value from all predictions 

p <-  ggplot(dat_expr, aes(x = min_reads, y =recall, color = expression)) + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + ggtitle("Exons with >0 mapped junction reads")
ggsave(paste0(OUTDIR, "rec_reads_expression.png"), p, device = "png")

p <- ggplot(dat_expr, aes(x = recall, y = precision, color = expression))  + geom_path(size = 1.5) + theme_bw() + ylim(c(0, 1)) + xlim(c(0,1)) + ggtitle("Exons with >0 mapped junction reads") + ylab("precision (all mapped exons)")
ggsave(paste0(OUTDIR, "PR_expression.png"), p, device = "png")


