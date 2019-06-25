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
    geom_col(position = position_dodge2( preserve = "total", padding=0.1), size=0.8) + 
    facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
    theme_bw() +
    theme(legend.position="bottom") +
    scale_y_log10() +
    guides(fill = guide_legend(nrow=2,byrow=TRUE), color = guide_legend(nrow=1)) +
    theme(legend.box = "vertical") +
    scale_color_manual(values=c("gold", "grey"))+
    theme(strip.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
    
  
  ggsave(file.path(out_dir,  paste0(prefix, "_evaluation_SJ_barplot.pdf")), 
         p, width = 7.5, height = 7)
  
  ## compute percentage
  dat <- dat %>%
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(total = sum(count),
                  percentage = count/total) %>%
    dplyr::ungroup()

  p <- ggplot(dat, aes(x = measure, y = percentage, 
                       fill = parameter, color = mapper)) +
    geom_col(position = position_dodge2( preserve = "total", padding=0.1), size=0.8) +
    facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
    theme_bw() +
    theme(legend.position="bottom") +
    guides(fill = guide_legend(nrow=2,byrow=TRUE), color = guide_legend(nrow=1)) +
    theme(legend.box = "vertical") +
    scale_color_manual(values=c("gold", "grey"))+
    theme(strip.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14)) +
    scale_y_sqrt()
  ggsave(file.path(out_dir, 
                   paste0(prefix, "_evaluation_SJ_barplot_percent.pdf")), 
         p, width = 7.5, height = 7) 
   
  ## plot the number of mapped reads per BAM file
  p <- ggplot(dat[dat$measure == "TP" & dat$read == "first",], 
              aes(x = mapper, y= total, fill = parameter) ) +
    geom_col(position = position_dodge2( preserve = "total", padding=0.1)) +
    facet_grid(rows = vars(removed_annotation)) +
    theme_bw() +
    theme(legend.position="bottom") +
    geom_text(aes(label = total, y = total/10, vjust = 0.5, hjust = 0), 
              position = position_dodge(0.9), angle = 80) +
    theme(strip.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  ggsave(file.path(out_dir, paste0(prefix, "_BAM_unique_mapped_barplot.pdf")), 
         p, width = 7, height = 7) 
  
  #############
  ## Accurracy: (TP + TN) / (TP + TN + FP + FN)
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
    geom_text(aes(label = round(accuracy, digits = 3), y = 0.1, 
                  vjust = 0.5, hjust = 0), 
              position = position_dodge(0.9), angle = 90) +
    theme(strip.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  ggsave(file.path(out_dir, paste0(prefix, 
                                   "_evaluation_SJ_barplot_accuracy.pdf")), 
         p, width = 7, height = 7) 
  
  ############
  ## ROC plot: FPR  vs TPR (sensitivity)
  dat <- dat %>% 
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(FPR = count[measure=="FP"]/(count[measure=="TN"] + 
                                              count[measure=="FP"]),
                  TPR = count[measure=="TP"]/(count[measure=="TP"] + 
                                              count[measure=="FN"]))
  

  max_axis <- round(max(max(dat$FPR), 1-min(dat$TPR)), digits = 2) + 0.01
  p <- ggplot(dat[dat$measure == "TP",], aes(x = FPR, y = TPR, 
                                             color = parameter)) +
    geom_point(aes(shape = mapper), size = 3) + 
    facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
    theme_bw() +
    scale_x_continuous(limits = c(0, max_axis)) +
    scale_y_continuous(limits = c(1-max_axis, 1)) +
    coord_equal(ratio=1) +
    theme(strip.text = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
  
  ggsave(file.path(out_dir, paste0(prefix, "_evaluation_SJ_ROC.pdf")), 
         p, width = 7, height = 6) 
  
  
  ############
  ## Precision-Recall plot: Recall (TPR) vs Precision
  dat <- dat %>% 
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(Precision = count[measure=="TP"]/(count[measure=="TP"] + 
                                                count[measure=="FP"]))
  
    p <- ggplot(dat[dat$measure == "TP",], aes(x = TPR, y = Precision, 
                                               color = parameter)) +
      geom_point(aes(shape = mapper), size = 3) + 
      facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
      theme_bw() +
      scale_x_continuous(limits = c(0.95, 1)) +
      scale_y_continuous(limits = c(0.95, 1)) +
      coord_equal(ratio=1) +
      labs(x = "Recall") +
      theme(strip.text = element_text(size = 12), 
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
    

    ggsave(file.path(out_dir, paste0(prefix, "_evaluation_SJ_PR.pdf")), 
           p, width = 7, height = 6) 

  #############
  ## F1 score: 2*(Precision*Recall)/(Precision+Recall)
  dat <- dat %>%
    dplyr::group_by(removed_annotation, parameter, mapper, read) %>%
    dplyr::mutate(F1 = 2 * (Precision * TPR) / (Precision + TPR))
    
    p <- ggplot(dat[dat$measure == "TP",], aes(x = mapper, y = F1, 
                                               fill = parameter )) + 
      geom_col(position = position_dodge2( preserve = "total", padding=0.1)) +
      facet_grid(rows = vars(removed_annotation), cols = vars(read)) +
      theme_bw() +
      theme(legend.position="bottom") +
      # geom_text(aes(label = round(F1, digits = 3), y = 0., 
      #               vjust = 0.5, hjust = 0), 
      #           position = position_dodge(0.9), angle = 90) +
      labs(y = "F1 score") +
      theme(strip.text = element_text(size = 12), 
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14)) +
      coord_cartesian(ylim=c(0.96, 1))
    
    ggsave(file.path(out_dir, paste0(prefix, 
                                     "_evaluation_SJ_barplot_F1.pdf")), 
           p, width = 7, height = 7) 
    
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
