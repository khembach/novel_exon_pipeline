args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

print(GTF)
print(SJFILE)
print(BAM)
print(OUTFILE)
print(OVERHANGMIN)

suppressPackageStartupMessages({
    library(exondiscovery)
})

## Prepare annotation
anno <- prepare_annotation(GTF)

## Predict novel exons using the exondiscovery package
novel_exon_df <- find_novel_exons(sj_filename = SJFILE, 
                                  annotation = anno, 
                                  min_unique = 1, 
                                  bam = BAM,
                                  overhang_min = as.integer(OVERHANGMIN))

## Write results
write.table(novel_exon_df, file = OUTFILE, row.names = FALSE, quote = FALSE, sep = "\t")
