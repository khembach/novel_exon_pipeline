args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

print(GTF)
print(EXONPRED)
print(OUTFILE)

suppressPackageStartupMessages({
    library(data.table)
    library(exondiscovery)
    library(rtracklayer)
})


novel_exon_df <- fread(EXONPRED)
gtf_added <- extend_gtf(GTF, novel_exon_df)
export(object = gtf_added, con = OUTFILE)
