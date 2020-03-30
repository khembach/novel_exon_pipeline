### Plot precision recall curves for discerns paper.
## First, we evaluate the exons by class. Then we only evaluate the exons with at least 1 mapped junction read (all others are not expected to be predicted). 
## Last, we we separate the exons according to their location, type or expression.


## We compute the two dat_class_expr dataframe for DISCERNS and StringTie (using the best paramter sets that we want to compare)
# we combine the two dataframes and plot PR curves split by class

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(data.table)
})

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}


## ---------------------------------------------------------------------------


evaluate_prediction <-function(criteria, pred, truth){
  ## check if the predictions can be found among the truth and add a column "label" with TRUE, otherwise NA
  labeled <- left_join( pred %>% select(criteria),  truth %>% select(c(criteria, "ID")) %>% mutate(label = TRUE) ) %>% select(matched_ID = ID, label)
  labeled$label[ is.na(labeled$label) ] <- FALSE
  cbind(pred, labeled)
}

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

## given a table with labeled predicted exons, compute precision and recall
compute_rec <- function(novelExons, nr_removed){
  ## compute TP and FP for each min_read bin
  eval_dat <- novelExons %>% group_by(min_reads) %>%
    summarise(tp = sum(label)) %>% arrange(desc(min_reads)) %>%
    mutate(sum_tp = cumsum(tp)) 
  ## compute FN
  eval_dat <- eval_dat %>% mutate(fn = nr_removed - sum_tp) 
  ## REC
  eval_dat <- eval_dat %>% mutate(recall = sum_tp/(sum_tp + fn))
  eval_dat
}

## ---------------------------------------------------------------------------


REMOVED <- "../simulation/analysis/mapped_junction_count/removed_me_exon_unique_classified_outSJfilterOverhangMin6_junc_count.txt"
PREDICTION <- "../simulation/analysis/filtered_SJ/package/me_exon/novel_exons_outSJfilterOverhangMin6.txt"
OUTDIR <- "../simulation/analysis/exon_prediction_performance/PR/discerns_paper_figures/"


PREDICTION_stringtie <- "../simulation/analysis/stringtie/predictions/me_exon/minReadCoverage1_minIsoformAbundance0.05/novel_exons_outSJfilterOverhangMin6_stringtie.txt"





## ---------------------------------------------------------------------------
compute_pr_tab <- function(PREDICTION, REMOVED){

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


  ###### add class and true junc counts, true gene count and add the matched truth ID
  ## label all predictions with TRUE or FALSE depending if the prediction is correct or not
  ## We include lend for left_terminal and rstart for right_terminal evaluation, because they should be NA both in the prediction and the removed list

  cassette <- c("seqnames", "start", "end", "strand", "lend", "rstart")
  left_terminal <- c("seqnames",  "end", "strand", "lend", "rstart")
  right_terminal <- c("seqnames",  "start", "strand", "lend", "rstart")
  eval_criteria <- list(cassette = cassette, left_terminal = left_terminal, right_terminal = right_terminal)

  novelExons_split <- split(novelExons, novelExons$location)
  novelExons_split <- lapply(seq_along(eval_criteria), function(x) evaluate_prediction(eval_criteria[[x]], novelExons_split[[names(eval_criteria)[x]]], removed))
  novelExons <- unsplit(novelExons_split,  novelExons$location)


  ####### We only consider the exons with at least 1 junction read
  ## Seperate the exons according to their class
  ## We do not know from which class the predictions are. Thus, we compute local recall but not precision (requires number of FP). As precision, we use the global precision (precision for a given number of supporting reads)

  keep <- which(removed$location == "cassette" & 
    pmin(removed$count_lsj, removed$count_rsj) > 0)
  keep <- c(keep,
            which(removed$location == "left_terminal" & removed$count_rsj > 0))
  keep <- c(keep,
            which(removed$location == "right_terminal" & removed$count_lsj > 0))
  removed <- removed[keep,]
  # removed <- removed %>% filter(pmin(removed$count_lsj, removed$count_rsj, na.rm = TRUE) > 0)  ## all exons with at least 1 junction read

  dat_all <- compute_pr_rec(novelExons, nrow(removed))
  precision_all <- dat_all %>% select(min_reads, precision)

  ###### me, normal exons
  id_type <- split(removed$ID, removed$type)  ## all IDs that belong to a certain class
  novelExons <- novelExons %>% mutate(type=ifelse((end-start+1)<=27, "me", "exon"))

  ###### low and high expressed exons: simulated number of reads below or above median
  median_expression <- median(removed$count_reads)

  ### low and high expression, AND different classes
  id_class_expr <- split(removed$ID, removed$class)  ## all IDs that belong to a certain class
  id_class_expr <- lapply( id_class_expr, function(x)
  split(x, ifelse(removed[ match(x, removed$ID), "count_reads"] < median_expression, "low", "high") )
  )

  dat_class_expr <- list()
  for(cl in names(id_class_expr)){
    for(exp in names(id_class_expr[[cl]])){
      dat_class_expr[[cl]][[exp]] <- novelExons %>% filter(matched_ID %in% id_class_expr[[cl]][[exp]])
    }
  }
  dat_class_expr <- lapply(names(id_class_expr), 
                           function(x) lapply(names(dat_class_expr[[x]]), 
                              function(u) compute_rec(dat_class_expr[[x]][[u]], 
                                                      length(id_class_expr[[x]][[u]]) )))

  names(dat_class_expr) <- names(id_class_expr)
  for(cl in names(id_class_expr)){
    names(dat_class_expr[[cl]]) <- names(id_class_expr[[cl]])
  }
  dat_class_expr <- lapply(dat_class_expr, function(x) rbindlist(x, idcol="expr" )  )
  dat_class_expr <- rbindlist(dat_class_expr, idcol = "class")
  dat_class_expr <- left_join(dat_class_expr, precision_all) ## add the corresponding precision value from all predictions 
  dat_class_expr
}



## DISCERNS predictions
discerns <- compute_pr_tab(PREDICTION, REMOVED)

p <- ggplot(discerns, aes(x = recall, y = precision, color = class))  + geom_path(size = 1.5, aes(linetype = expr ) ) + theme_bw(base_size = 14) + geom_point()+ ylim(c(0, 1)) + xlim(c(0,1))+ ggtitle("Exons with >0 mapped junction reads")+ ylab("precision (all mapped exons)") + 
scale_linetype_manual(values=c("solid", "twodash"))
ggsave(paste0(OUTDIR, "PR_class_expr_discerns.pdf"), p, device = "pdf", height = 5, width = 6)


## StringTie predictions
stringtie <- compute_pr_tab(PREDICTION_stringtie, REMOVED)

p <- ggplot(stringtie, aes(x = recall, y = precision, color = class))  + 
geom_path(size = 1.5, aes(linetype = expr ) ) + 
theme_bw(base_size = 14) + 
geom_point()+ 
ylim(c(0, 1)) + 
xlim(c(0,1))+ 
ggtitle("Exons with >0 mapped junction reads")+ 
ylab("Precision (all mapped exons)") + 
scale_linetype_manual(values=c("solid", "twodash"))

ggsave(paste0(OUTDIR, "PR_class_expr_stringtie.pdf"), p, device = "pdf", height = 5, width = 6)


###### 
# Comparison
discerns$method <- "DISCERNS"
stringtie$method <- "StringTie"

dat <- rbind(discerns, stringtie)
dat$class <- factor(dat$class, levels = c("easy", "medium", "complicated"))
dat$expr <- factor(dat$expr, levels = c("low", "high"))

dat$group <- paste0(dat$method, dat$expr)

p <- ggplot(dat, aes(x = recall, y = precision, color = method)) + 
  geom_path(aes(group = group ), size = 0.8) + 
  theme_bw(base_size = 14) + 
  geom_point(aes(shape = expr), size = 2) + 
  ylim(c(0, 1)) + 
  xlim(c(0, 1)) + 
  xlab("Recall") +
  ylab("Precision") + 
  facet_wrap(~class) + 
  coord_equal(ratio=1) + 
  theme(legend.position="bottom", legend.box = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values=c(16, 15)) +
  labs(shape = "expression") +
  scale_color_manual(values = c("#332288", "#F17E06"))

# ggsave(paste0(OUTDIR, "PR_class_expr_comparison.pdf"), p, device = "pdf", height = 4, width = 7.5)
ggsave(paste0(OUTDIR, "PR_class_expr_comparison.png"), p, device = "png",
height = 3.5, width = 6.3, dpi = 300)


