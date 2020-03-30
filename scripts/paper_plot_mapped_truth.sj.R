## This is R script generates the plots for the mapping comparison (STAR vs. hisat2) in the discerns paper.


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
eval_sj <- function(in_files, in_path) {
  ###############
  ## Read the offset files
  eval_tabs <- lapply(in_files, function(x) read_eval_measures(x, in_path))
  dat <- do.call("rbind", eval_tabs)
  dat$mapper <- as.factor(dat$mapper)
  dat$removed_annotation <- as.factor(dat$removed_annotation)
  dat$parameter <- as.factor(dat$parameter)

  ############
  ## ROC plot: FPR  vs TPR (sensitivity)
  dat <- dat %>% 
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(FPR = count[measure=="FP"]/(count[measure=="TN"] + 
                                              count[measure=="FP"]),
                  TPR = count[measure=="TP"]/(count[measure=="TP"] + 
                                              count[measure=="FN"]))
  ############
  ## Precision-Recall plot: Recall (TPR) vs Precision
  dat <- dat %>% 
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(Precision = count[measure=="TP"]/(count[measure=="TP"] + 
                                                count[measure=="FP"]))
  
  return(dat)
}



plot_PR_curve <- function(dat, out_dir, prefix){
  ## TODO: legend at bottom, bigger text, bigger and transparent symbols
  ## combine plot for all and SJ reads
    p <- ggplot(dat[dat$measure == "TP",], aes(x = TPR, y = Precision, 
                                               color = parameter)) +
      geom_point(aes(shape = mapper), size = 4, alpha = 0.7) + 
      facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
      theme_bw(base_size = 14) +
      scale_x_continuous(limits = c(0.95, 1)) +
      scale_y_continuous(limits = c(0.95, 1)) +
      coord_equal(ratio=1) +
      labs(x = "Recall") +
      theme(strip.text = element_text(size = 14), 
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            legend.position="bottom", legend.box = "vertical") +
      guides(shape = guide_legend(order = 1),
             color = guide_legend(order = 2, ncol = 2))

    ggsave(file.path(out_dir, paste0(prefix, "_evaluation_SJ_PR.pdf")), 
           p, width = 7, height = 7) 
}


plot_PR_curve_combined <- function(dat1, dat2, out_dir, prefix){
  dat <- rbind(dat1 %>% mutate(reads = "all reads"), dat2 %>% mutate(reads = "novel SJ reads"))


## We sum up the measures for first and second reads and recompute Precision and Recall
  dat <- dat %>% 
    group_by(measure, mapper, removed_annotation, parameter, reads) %>%
    summarise(count = sum(count))

  dat <- dat %>% 
    dplyr::group_by(removed_annotation, parameter, mapper, reads) %>%
    dplyr::mutate(FPR = count[measure=="FP"]/(count[measure=="TN"] + 
                                              count[measure=="FP"]),
                  TPR = count[measure=="TP"]/(count[measure=="TP"] + 
                                              count[measure=="FN"]),
                  Precision = count[measure=="TP"]/(count[measure=="TP"] + 
                                                count[measure=="FP"]),
                  F1 = 2 * (Precision * TPR) / (Precision + TPR))
  levels(dat$removed_annotation) <- gsub("_", " + ", levels(dat$removed_annotation))


  pr <- ggplot(dat[dat$measure == "TP",], aes(x = TPR, y = Precision, 
                                             color = parameter)) +
    geom_point(aes(shape = mapper), size = 4, alpha = 0.7) + 
    facet_grid(rows = vars(removed_annotation), cols = vars(reads), labeller = label_wrap_gen(width = 10)) +
    theme_bw(base_size = 14) +
    scale_x_continuous(limits = c(0.95, 1)) +
    scale_y_continuous(limits = c(0.95, 1)) +
    coord_equal(ratio=1) +
    labs(x = "Recall") +
    theme(strip.text = element_text(size = 14), 
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          legend.position = "none") +
    scale_color_locuszoom()
    # scale_color_uchicago()
    # scale_color_d3()
    #       legend.position="bottom", legend.box = "vertical") +
    # guides(shape = guide_legend(order = 1),
    #        color = guide_legend(order = 2, ncol = 2))

  ggsave(file.path(out_dir, paste0(prefix, "_evaluation_SJ_PR.pdf")), 
         pr, width = 5, height = 6) 

## F1 score plot
  dat <- dat %>% mutate(combination = paste0(mapper, "_", parameter))
  p <- ggplot(dat[dat$measure == "TP",], aes(x = combination, y = F1, 
                                             color = parameter, shape = mapper )) + 
    geom_segment(aes(x = combination, xend = combination, y = 0.95, yend = F1)) +
    geom_point(size=4, alpha=0.6) +
    theme_bw(base_size = 14) +
    facet_grid(rows = vars(removed_annotation), cols = vars(reads)) +
    labs(y = "F1 score", x = "Mapping tool and parameter") +
    theme(strip.text = element_text(size = 14), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.title = element_text(size = 16),
          legend.position = "none") +
    scale_y_continuous(limits = c(0.95, 1)) + 
    scale_color_locuszoom() +
    coord_flip()
  
  ggsave(file.path(out_dir, paste0(prefix, 
                                   "_evaluation_SJ_lollipop_F1.pdf")), 
         p, width = 5, height = 6) 

  ## plot only the legend
  p <- ggplot(dat, aes(x = TPR, y = Precision, color = parameter, 
                       shape = mapper)) + 
    geom_point(size = 4) + 
    theme_bw(base_size = 20) +
    theme(legend.position="bottom", legend.box = "horizontal") +
    guides(shape = guide_legend(order = 1),
           color = guide_legend(order = 2, ncol = 2)) +
    scale_color_locuszoom()


  legend <- cowplot::get_legend(p)
  pdf(file.path(out_dir, paste0(prefix,"_legend.pdf")), width = 12, height = 1.5)
  grid.newpage()
  grid.draw(legend)
  dev.off()
 

}


##----------------------------------------------------------------------------
## 


library(ggplot2)
library(dplyr)
library("ggsci")
library(cowplot)
library(grid)

## Parameters
# INDIR <- snakemake@input[["indir"]]
# OUTDIR <- snakemake@params[["outdir"]]

INDIR <- "../simulation/mapped_truth"
OUTDIR <- "../simulation/analysis/mapped_sj_eval/discerns_paper_figures"

in_path <- INDIR
out_dir <- OUTDIR
prefix <- "all_reads"


#############
## Plotting
print("all reads")
in_files <- list.files(INDIR, pattern = "_evaluation_SJ_all.txt", 
                           recursive = TRUE)
dat_all <- eval_sj(in_files, INDIR)
plot_PR_curve(dat_all, OUTDIR, prefix = "all_reads")

print("removed exons reads")
in_files_removed <- list.files(INDIR, 
                           pattern = "_evaluation_SJ_overl_removed_exons.txt", 
                           recursive = TRUE)
dat_removed <- eval_sj(in_files_removed, INDIR)
plot_PR_curve(dat_removed, OUTDIR, prefix = "reads_removed_exons")


## both curves in one pdf
plot_PR_curve_combined(dat_all, dat_removed, OUTDIR, prefix = "all_and_removed_exons")


