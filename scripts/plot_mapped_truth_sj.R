## This script computes the accuracy per BAM file and plots the performance of
## aligned and paramter combination as a point on a ROC plot.

#' Split a file path into folder names
#' 
#' Split a file path into a character vector of components
#'
#' @param path file path
#'
#' @return character vector with folder names
#' @export
#'
#' @examples
split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

#' Read evaluation measures
#' 
#' Read the evaluation measures from file and return a data.frame with the 
#' measures, the mapper, the removed annotation and the specific parameter(s).
#'
#' @param file path to the input file
#' @param path relative path to the file
#'
#' @return data.frame with the evaluation measures, the used mapper, the removed
#' annotation and the specific parameter
#' @export
#'
#' @examples
read_eval_measures <- function(file, path){
  ## split the folder names
  dir_names <- rev( split_path(file) )
  dir_names <- dir_names[-length(dir_names)]
  
  df <- read.table(file.path(path, file), header = TRUE, sep = "\t")
  df$mapper <- dir_names[1]
  df$removed_annotation <- dir_names[2]
  df$parameter <- if(length(dir_names)==3) dir_names[3] else "default"
  df
}


#' Plots from SJ evaluation results
#' 
#' Create plots with number of reads per BAM file, read count and percentage per
#' measure, accuracy and ROC plot.
#'
#' @param in_files vector with paths to the SJ evaluation files (relative to in_path)
#' @param in_path path to the base directory of the SJ evaluation files 
#' @param outdir path to the output directory in which the plots will be saved
#' @param prefix output file prefix
#'
#' @return
#' @export
#'
#' @examples
eval_sj <- function(in_files, in_path, out_dir, prefix) {
  ###############
  ## Read the offset files
  eval_tabs <- lapply(in_files, function(x) read_eval_measures(x, in_path))
  dat <- do.call("rbind", eval_tabs)
  dat$mapper <- as.factor(dat$mapper)
  dat$removed_annotation <- as.factor(dat$removed_annotation)
  dat$parameter <- as.factor(dat$parameter)

  p <- ggplot(dat, aes(x = measure, y = count, 
                       fill = parameter, color = mapper)) +
    geom_col(position = position_dodge2( preserve = "total", padding=0.1)) + 
    facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
    theme_bw() +
    theme(legend.position="bottom") +
    scale_y_log10() +
    guides(fill = guide_legend(nrow=2,byrow=TRUE), color = guide_legend(nrow=1)) +
    theme(legend.box = "vertical") 
  
  ggsave(file.path(out_dir,  paste0(prefix, "_evaluation_SJ_barplot.pdf")), 
         p, width = 7, height = 7)
  
  ## compute percentage
  dat <- dat %>%
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(total = sum(count),
                  percentage = count/total) %>%
    dplyr::ungroup()

  p <- ggplot(dat, aes(x = measure, y = percentage, 
                       fill = parameter, color = mapper)) +
    geom_col(position = position_dodge2( preserve = "total", padding=0.1)) +
    facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
    theme_bw() +
    theme(legend.position="bottom") +
    guides(fill = guide_legend(nrow=2,byrow=TRUE), color = guide_legend(nrow=1)) +
    theme(legend.box = "vertical") 
  ggsave(file.path(out_dir, 
                   paste0(prefix, "_evaluation_SJ_barplot_percent.pdf")), 
         p, width = 7, height = 7) 
   
  ## plot the number of mapped reads per BAM file
  p <- ggplot(dat[dat$measure == "TP" & dat$read == "first",], 
              aes(x = mapper, y= total, fill = parameter) ) +
    geom_col(position = position_dodge2( preserve = "total", padding=0.1)) +
    facet_grid(rows = vars(removed_annotation)) +
    theme_bw() +
    theme(legend.position="bottom") +
    geom_text(aes(label = total, y = total/10, vjust = 0.5, hjust = 0), 
              position = position_dodge(0.9), angle = 80)
  
  ggsave(file.path(out_dir, paste0(prefix, "_BAM_unique_mapped_barplot.pdf")), 
         p, width = 7, height = 7) 
  
  #############
  ## Accurracy
  dat <- dat %>%
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(accuracy = (count[measure=="TP"] + count[measure=="TN"])/
                    (count[measure=="TP"] + count[measure=="TN"] + 
                     count[measure=="FP"] + count[measure=="FN"]))
  
  
  p <- ggplot(dat[dat$measure == "TP",], aes(x = mapper, y = accuracy, 
                                             fill = parameter )) + 
    geom_col(position = position_dodge2( preserve = "total", padding=0.1)) +
    facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
    theme_bw() +
    theme(legend.position="bottom") +
    geom_text(aes(label = round(accuracy, digits = 4), y = 0.1, 
                  vjust = 0.5, hjust = 0), 
              position = position_dodge(0.9), angle = 85)
  
  ggsave(file.path(out_dir, paste0(prefix, 
                                   "_evaluation_SJ_barplot_accuracy.pdf")), 
         p, width = 7, height = 7) 
  
  ############
  ## ROC plot
  dat <- dat %>% 
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(FPR = count[measure=="TP"]/(count[measure=="TP"] + 
                                              count[measure=="FP"]),
                  TPR = count[measure=="TP"]/(count[measure=="TP"] + 
                                              count[measure=="FN"]))
  
  p <- ggplot(dat[dat$measure == "TP",], aes(x = FPR, y = TPR, 
                                             color = parameter)) +
    geom_point(aes(shape = mapper)) + 
    facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
    theme_bw() +
    scale_x_continuous(limits = c(0.95, 1)) +
    scale_y_continuous(limits = c(0.95, 1)) +
    coord_equal(ratio=1)
    
  ggsave(file.path(out_dir, paste0(prefix, "_evaluation_SJ_barplot_ROC.pdf")), 
         p, width = 7, height = 6) 
}


##----------------------------------------------------------------------------
## 

library(ggplot2)
library(dplyr)

## Parameters
INDIR <- snakemake@input[["indir"]]
OUTDIR <- snakemake@params[["outdir"]]

# INDIR <- "../simulation/mapped_truth"
# OUTDIR <- "../simulation/analysis/mapped_sj_eval/"

#############
## Plotting
print("all reads")
in_files <- list.files(INDIR, pattern = "_evaluation_SJ_all.txt", 
                           recursive = TRUE)
eval_sj(in_files, INDIR, OUTDIR, prefix = "all_reads")

print("removed exons reads")
in_files <- list.files(INDIR, 
                           pattern = "_evaluation_SJ_overl_removed_exons.txt", 
                           recursive = TRUE)
eval_sj(in_files, INDIR, OUTDIR, prefix = "reads_removed_exons")
