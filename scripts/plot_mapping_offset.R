### This script plots the mapping offset of different samples.
### and the sum for each of the different mappers, parameters, gtf files

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

#' Read offsets from file
#' 
#' Read the offset table from file and return a data.frame with the offset counts, 
#' the mapper, the removed annotation and the specific parameter(s).
#'
#' @param file path of the offset count table, relative to the mapped_truth directory
#'
#' @return data.frame with the offset counts, the used mapper, the removed annotation
#' and the specific parameter
#' @export
#'
#' @examples
read_offsets <- function(file){
  ## split the folder names
  dir_names <- rev( split_path(file) )
  dir_names <- dir_names[-length(dir_names)]
  
  offsets <- read.table(file.path(OFFSET_DIR, file), header = TRUE, sep = "\t")
  offsets$mapper <- dir_names[1]
  offsets$removed_annotation <- dir_names[2]
  offsets$parameter <- if(length(dir_names)==3) dir_names[3] else "default"
  
  offsets
}


#' Offset distribution plots
#' 
#' Create plots to show the distribution of offsets in the hisat2 and star alignments.
#'
#' @param offset_files vector with paths to the offset tables
#' @param outdir path to the output directory in which the plots will be saved
#'
#' @return
#' @export
#'
#' @examples
plot_offset_distribution <- function(offset_files, outdir, prefix) {
  ###############
  ## Read the offset files
  offset_tabs <- lapply(offset_files, read_offsets)
  dat <- do.call("rbind", offset_tabs)
  dat$mapper <- as.factor(dat$mapper)
  dat$removed_annotation <- as.factor(dat$removed_annotation)
  dat$parameter <- as.factor(dat$parameter)
  
  ###########
  ## Plot the offset distributions
  p <- ggplot(dat, aes(x=offset, y=count, color = parameter)) +
    geom_line(alpha=0.5) + 
    facet_grid(rows = vars(removed_annotation), cols = vars(mapper), scales="fixed") +
    theme_bw() +
    scale_y_log10() +
    theme(legend.position="bottom")
  ggsave(file.path(OUTDIR, paste0(prefix, "_mapped_offset_comparison.png")), p, width=7, height = 6)
  
  ## compute the sum of wrong reads (offset >=101), the shifted reads (offset between
  ## 1 and 100) and the correct reads (offset == 0)
  dat_sums <- dat %>% 
    dplyr::group_by(removed_annotation, parameter, mapper) %>% 
    dplyr::summarize(correct = count[offset == 0], 
                     shifted = sum(count[offset>0 & offset<101]), 
                     wrong = count[offset == 101])
  ## transform to long format
  dat_sums <- dat_sums %>% 
    gather(key = read_mapping, value = count, correct, shifted, wrong, factor_key = TRUE)
  ## compute percentage
  dat_sums <- dat_sums %>%
    dplyr::group_by(removed_annotation, parameter, mapper) %>%
    dplyr::mutate(total = sum(count),
                  percentage = count/total)
  
  ## write the count tabel to file
  write.table(dat_sums, file.path(OUTDIR,  paste0(prefix, "_read_offset_table.txt")), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ###########
  ## Barplot
  p <- ggplot(dat_sums, aes(x = read_mapping, y = count, fill = parameter)) +
    geom_col(position = "dodge") +
    facet_grid(rows = vars(removed_annotation), cols = vars(mapper)) +
    scale_y_log10() + 
    theme_bw() +
    theme(legend.position="bottom")
  ggsave(file.path(OUTDIR,  paste0(prefix, "_mapped_offset_count_barplot.png")), p, width = 7, height = 7)
  
  p <- ggplot(dat_sums, aes(x = read_mapping, y = percentage, fill = parameter)) +
    geom_col(position = "dodge") +
    facet_grid(rows = vars(removed_annotation), cols = vars(mapper)) +
    theme_bw() +
    theme(legend.position="bottom")
  ggsave(file.path(OUTDIR,  paste0(prefix, "_mapped_offset_perc_barplot.png")), p, width = 7, height = 7)
  
}


###############
## Parameters
OFFSET_DIR <- snakemake@input[["offset_dir"]]
OUTDIR <- snakemkae@output[["outdir"]]

# OFFSET_DIR <- "../simulation/mapped_truth"
# OUTDIR <- "../simulation/analysis/mapped_offset/"

library(ggplot2)
library(dplyr)
library(tidyr)

#############
## Plotting
offset_files <- list.files(OFFSET_DIR, pattern = "mapped_truth.txt", recursive = TRUE)
plot_offset_distribution(offset_files, OUTDIR, prefix = "all_reads")

offset_files <- list.files(OFFSET_DIR, pattern = "mapped_truth_removed_exons.txt", recursive = TRUE)
plot_offset_distribution(offset_files, OUTDIR, prefix = "reads_removed_exons")
