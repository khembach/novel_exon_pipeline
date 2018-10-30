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
#' @param prefix prefix of the ouput files
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
  ggsave(file.path(outdir, paste0(prefix, "_mapped_offset_comparison.pdf")), p, width=7, height = 6)
  
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
  
  ## write the count table to file
  write.table(dat_sums, file.path(outdir,  paste0(prefix, "_read_offset_table.txt")), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ###########
  ## Barplot
  p <- ggplot(dat_sums, aes(x = read_mapping, y = count, fill = parameter)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single", padding=0.05)) +
    facet_grid(rows = vars(removed_annotation), cols = vars(mapper)) +
    scale_y_log10() + 
    theme_bw() +
    theme(legend.position="bottom")
  ggsave(file.path(outdir,  paste0(prefix, "_mapped_offset_count_barplot.pdf")), p, width = 7, height = 7)
  
  p <- ggplot(dat_sums, aes(x = read_mapping, y = percentage, fill = parameter)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single", padding=0.05)) +
    facet_grid(rows = vars(removed_annotation), cols = vars(mapper)) +
    theme_bw() +
    theme(legend.position="bottom")
  ggsave(file.path(outdir,  paste0(prefix, "_mapped_offset_perc_barplot.pdf")), p, width = 7, height = 7)
  
}



#' Plots of offset corrected by soft clipped bases
#' 
#' Create plots to show the relationship between the offset and the number of soft
#' clipped bases per read.
#'
#' @param offset_files vector with paths to the offset_soft_clipped tables
#' @param outdir path to the output directory in which the plots will be saved
#' @param prefix prefix of the ouput files
#'
#' @return
#' @export
#'
#' @examples
plot_offset_sc <- function(offset_files, outdir, prefix) {
  ###############
  ## Read the offset files
  offset_tabs <- lapply(offset_files, read_offsets)
  dat <- do.call("rbind", offset_tabs)
  dat$mapper <- as.factor(dat$mapper)
  dat$removed_annotation <- as.factor(dat$removed_annotation)
  dat$parameter <- as.factor(dat$parameter)

  ###########
  ## Plot the offset and the corresponding number of soft clipped bases
  ## We only plot the reads with offset >0 and <=101
  dat_part <- dat[dat$parameter=="default" & dat$offset<=101 & dat$offset>0, ]
  
  ## binhex plot
  p <- ggplot(dat_part, aes(x=offset, y=soft_clipped)) +
    geom_hex(bins=50) + 
    facet_grid(rows = vars(removed_annotation), cols = vars(mapper), scales="fixed") +
    theme_bw() +
    theme(legend.position="bottom") +
    labs(y = "number of soft clipped bases", title="All reads with 0 < offset <= 101; default parameters")
  ggsave(file.path(outdir, paste0(prefix, "_offset_soft_clipped_default.pdf")), p, width=5, height = 6)

  dat_part <- dat[dat$mapper!="hisat2" & dat$offset<=101 & dat$offset>0, ]
  p <- ggplot(dat_part, aes(x=offset, y=soft_clipped)) +
    geom_hex(bins=50) + 
    facet_grid(rows = vars(removed_annotation), cols = vars(parameter), scales="fixed") +
    theme_bw() +
    theme(legend.position="bottom") +
    labs(y = "number of soft clipped bases", title="All reads with 0 < offset <= 101")
  ggsave(file.path(outdir, paste0(prefix, "_offset_soft_clipped_star_params.pdf")), p, width=10, height = 6)
  
  ## most offsets are 0 if we subtract the # of soft clipped nts
  # dat_part <- dat[dat$offset<=101 & dat$offset>0, ]
  # p <- ggplot(dat_part, aes( y = abs(offset - soft_clipped), fill = parameter)) +
  #   geom_boxplot() +
  #   facet_grid(rows = vars(removed_annotation), cols = vars(mapper), space = "free_x") +
  #   theme_bw() +
  #   theme(legend.position="bottom") +
  #   theme(axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())
    

  ## We subtract the soft clipped bases from the offset, sum the result and divide
  ## by the number of reads
  sum_corrected_offset <- dat %>%
    dplyr::group_by(mapper, parameter, removed_annotation) %>%
    dplyr::summarize(offset_no_sc = sum(as.numeric(abs(offset - soft_clipped)))/n())
  
  write.table(sum_corrected_offset, file.path(outdir,  paste0(prefix, "_offset_without_sc.txt")), quote = FALSE, sep = "\t", row.names = FALSE)
  
  p <- ggplot(sum_corrected_offset, aes(x=mapper, y=offset_no_sc, fill = parameter)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = "single", padding=0.05)) +
    facet_grid(rows = vars(removed_annotation), scales = "free_x") +
    theme_bw() +
    theme(legend.position="bottom") +
    scale_y_log10() +
    labs(y = "sum(offset - soft clipped) / N_reads")
  ggsave(file.path(outdir, paste0(prefix, "_offset_without_sc_barplot.pdf")), p, width=7, height = 6)
  
}


###############
## Parameters
OFFSET_DIR <- snakemake@input[["offset_dir"]]
OUTDIR <- snakemake@params[["outdir"]]

# OFFSET_DIR <- "../simulation/mapped_truth"
# OUTDIR <- "../simulation/analysis/mapped_offset/"

library(ggplot2)
library(dplyr)
library(tidyr)

#############
## Plotting
print("all reads")
offset_files <- list.files(OFFSET_DIR, pattern = "offset_counts.txt", recursive = TRUE)
plot_offset_distribution(offset_files, OUTDIR, prefix = "all_reads")
print("all reads, sc correction")
offset_files <-  list.files(OFFSET_DIR, pattern = "offset_soft_clipped.txt", recursive = TRUE)
plot_offset_sc(offset_files, OUTDIR, prefix = "all_reads")

print("removed exons reads")
offset_files <- list.files(OFFSET_DIR, pattern = "offset_counts_removed_exons.txt", recursive = TRUE)
plot_offset_distribution(offset_files, OUTDIR, prefix = "reads_removed_exons")
print("removed exons reads, sc correction")
offset_files <-  list.files(OFFSET_DIR, pattern = "offset_soft_clipped_removed_exons.txt", recursive = TRUE)
plot_offset_sc(offset_files, OUTDIR, prefix = "reads_removed_exons")


