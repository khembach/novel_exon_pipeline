library(dplyr)
library(ggplot2)


REMOVED <- "simulation/reduced_GTF/removed_microexons_exons_unique_classified.txt"
NOVEL <- "simulation/analysis/filtered_SJ/me_exon/novel_exons_outSJfilterOverhangMin6.txt"

params <- c("default", "outSJfilterOverhangMin9", "outSJfilterOverhangMin6", "outSJfilterCountTotalMin3", "scoreGenomicLengthLog2scale0", "alignSJoverhangMin3")
PRED_PATH <- "simulation/analysis/filtered_SJ/me_exon/novel_exons_"


OUTDIR <- "simulation/analysis/exon_prediction_performance/outSJfilterDistToOtherSJmin0/"

##  compare prediction to true removed exons
removed <- read.table(REMOVED, header=TRUE) ## 261 removed exons
## the unique removed exons, without the splice junctions
removed_unique <- removed[,-which(colnames(removed) %in% c("rstart", "lend"))]
removed_unique <- removed_unique[ rownames(unique(removed_unique[, c("seqnames", "start", "end", "strand")])),]

# table(is.na(removed$lend) | is.na(removed$rstart)) # 81 internal and 43 terminal me
# removed_internal <- removed[which(!is.na(removed$lend) & !is.na(removed$rstart)),]

removed_per_class <- split(removed, removed$class)

## performance results
res <- data.frame(param = rep(NA, 2*length(params)), annotation = rep(NA, 2*length(params)), TP = rep(NA, 2*length(params)), FP = rep(NA, 2*length(params)), FN = rep(NA, 2*length(params)), TPR = rep(NA, 2*length(params)), PPV = rep(NA, 2*length(params)))

ind <- 1
## for each mapping parameter:
for (p in params){
	print(p)
	pred <- read.table(paste0(PRED_PATH, p, ".txt"),header = TRUE)

	pred_unique <- pred[,-which(colnames(pred) %in% c("rstart", "lend"))]
	pred_unique <- pred_unique[ rownames(unique(pred_unique[, c("seqnames", "start", "end", "strand")])),]

	## compute TP predictions --> all four coordinates are correct (lend, start, end, rstart)
	tp <- dplyr::semi_join(removed, pred, by = c("seqnames", "start", "end", "strand", "lend", "rstart")) 
	tp_unique <- dplyr::semi_join(removed_unique, pred_unique, by = c("seqnames", "start", "end", "strand"))

	## compute FP predictions --> the four coordinates are not correct
	fp <- dplyr::anti_join(pred, removed, by = c("seqnames", "start", "end", "strand", "lend", "rstart"))
	fp_unique <- dplyr::anti_join(pred_unique, removed_unique, by = c("seqnames", "start", "end", "strand"))

	## compute FN predictions --> removed exons that were not correctly predicted
	fn <- dplyr::anti_join(removed, pred, by = c("seqnames", "start", "end", "strand", "lend", "rstart"))
	fn_unique <- dplyr::anti_join(removed_unique, pred_unique, by = c("seqnames", "start", "end", "strand"))


	tpr <- round( nrow(tp) / (nrow(tp) + nrow(fn) ), 3)
	tpr_unique <- round( nrow(tp_unique) / (nrow(tp_unique) + nrow(fn_unique) ), 3)

	ppv <- round( nrow(tp) / (nrow(tp) + nrow(fp) ), 3)
	ppv_unique <- round( nrow(tp_unique) / (nrow(tp_unique) + nrow(fp_unique) ), 3)

	fdr <- round( nrow(fp) / (nrow(tp) + nrow(fp) ), 3)
	fdr_unique <- round( nrow(fp_unique) / (nrow(tp_unique) + nrow(fp_unique) ), 3)

	## comput the values seperately for the easy, medium and complicated exons
	## easy = unique exons

	res[ind,] <- c(p, "all", nrow(tp), nrow(fp), nrow(fn), tpr, ppv)
	res[ind + 1,] <- c(p, "unique", nrow(tp_unique), nrow(fp_unique), nrow(fn_unique), tpr_unique, ppv_unique)

	# print("TP\tFP\tFN")
	# print(cbind(nrow(tp), nrow(fp), nrow(fn) ))
	# print(cbind(nrow(tp_unique), nrow(fp_unique), nrow(fn_unique) ))

	ind <- ind + 2
}



## df with results for the different exon classes
res_class <- data.frame(param = rep(NA, 3*length(params)), class = rep(NA, 3*length(params)), TP = rep(NA, 3*length(params)), FN = rep(NA, 3*length(params)), TPR = rep(NA, 3*length(params)) )

class <- unique(factor(removed$class))

ind <- 1

## for each mapping parameter:
for (p in params){
	print(p)
	pred <- read.table(paste0(PRED_PATH, p, ".txt"),header = TRUE)

	# pred_unique <- pred[,-which(colnames(pred) %in% c("rstart", "lend"))]
	# pred_unique <- pred_unique[ rownames(unique(pred_unique[, c("seqnames", "start", "end", "strand")])),]

	## comput the values seperately for the easy, medium and complicated exons
	## easy = unique exons
	for (cl in class){

		r <- removed_per_class[[cl]]

		## compute TP predictions --> all four coordinates are correct (lend, start, end, rstart)
		tp <- dplyr::semi_join(r, pred, by = c("seqnames", "start", "end", "strand", "lend", "rstart")) 

		## compute FN predictions --> removed exons that were not correctly predicted
		fn <- dplyr::anti_join(r, pred, by = c("seqnames", "start", "end", "strand", "lend", "rstart"))

		tpr <- round( nrow(tp) / (nrow(tp) + nrow(fn) ), 3)

		res_class[ind,] <- c(p, cl, nrow(tp), nrow(fn), tpr)

		# print("TP\tFP\tFN")
		# print(cbind(nrow(tp), nrow(fp), nrow(fn) ))
		# print(cbind(nrow(tp_unique), nrow(fp_unique), nrow(fn_unique) ))

		ind <- ind + 1
	}
}
res_class$TP <- as.numeric(res_class$TP)
res_class$FN <- as.numeric(res_class$FN)
res_class$TPR <- as.numeric(res_class$TPR)



### Write the results to files
write.table(res, file.path(OUTDIR, "performance_tpr_ppv.txt"), sep = "\t", row.names = FALSE, quote=FALSE )
write.table(res_class, file.path(OUTDIR, "performance_exon_class_tpr.txt"), sep = "\t", row.names = FALSE, quote=FALSE )

### plot the TPR for the different classes and parameters
g <- ggplot(data = res_class, aes(x = param, y = TPR, color = class)) + geom_point(size = 3) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1) )
g
ggsave(file.path(OUTDIR, "exon_class_tpr.pdf"))





## plot exon level ROC curve base on simulated number of reads per exon
## maybve plot ROC curve based on single sj, take number of reads with specific sj


##  compare prediction to true removed exons
removed <- read.table(REMOVED, header=TRUE) ## 261 removed exons
# table(is.na(removed$lend) | is.na(removed$rstart)) # 81 internal and 43 terminal me
removed_internal <- removed[which(!is.na(removed$lend) & !is.na(removed$rstart)),]

novelExons <- read.table(NOVEL, header=TRUE) ## 170



removed[,c("seqnames", "lend", "start", "end", "rstart", "strand")]
library(dplyr)
dplyr::setdiff(novelExons, removed[,c("seqnames", "lend", "start", "end", "rstart", "strand")]) ## FP predictions: 55

not_found <- dplyr::setdiff(removed[,c("seqnames", "lend", "start", "end", "rstart", "strand")], novelExons) ## all predictions are 
# dplyr::intersect(novelExons, removed[,c("seqnames", "lend", "start", "end", "rstart")])

predicted <- dplyr::semi_join(removed, novelExons) ## most are unique exons
missed <- dplyr::anti_join(removed, novelExons) ## most missed are overlapping or lowly expressed



dplyr::setdiff( novelExons[,c("seqnames", "lend", "start", "strand")], removed[,c("seqnames", "lend", "start", "strand")] )  ## FP predictions of first junction: 33

dplyr::setdiff( novelExons[,c("seqnames", "end", "rstart", "strand")], removed[,c("seqnames", "end", "rstart", "strand")] ) ## FP predictions of second junction: 35

dplyr::setdiff( novelExons[,c("seqnames",  "start", "end", "strand")], removed[,c("seqnames", "start", "end", "strand")] )  ## FP predictions of nvoel exon: 55
