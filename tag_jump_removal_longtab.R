#!/usr/bin/env Rscript

## Simplified Tag-jump removal script
## Inputs:
##   --seqtab : sequence table in long format (SampleID, SeqID, Abundance)
##   --uc     : UC parquet file from dereplication
##   --f      : f parameter for UNCROSS (default 0.01)
##   --p      : p parameter for UNCROSS (default 1)
##   --outdir : output directory


suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-s", "--seqtab"), action="store", type='character', help="Sequence table (long format)"),
  make_option(c("-u", "--uc"), action="store", type='character', help="UC parquet file from dereplication"),
  make_option(c("-f", "--uncross_f"), action="store", type='numeric', default=0.01, help="f parameter of UNCROSS"),
  make_option(c("-p", "--uncross_p"), action="store", type='numeric', default=1, help="p parameter of UNCROSS"),
  make_option(c("-o", "--outdir"), action="store", type='character', default=".", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

SEQ_TAB_FILE <- opt$seqtab
UC_FILE      <- opt$uc
F_PARAM      <- opt$uncross_f
P_PARAM      <- opt$uncross_p
OUTPUT_DIR   <- opt$outdir

cat("Sequence table:", SEQ_TAB_FILE, "\n")
cat("UC file:", UC_FILE, "\n")
cat("f parameter:", F_PARAM, "\n")
cat("p parameter:", P_PARAM, "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# Load sequence table
cat("Loading sequence table...\n")
SEQTAB <- fread(SEQ_TAB_FILE, header=FALSE)
setnames(SEQTAB, c("SampleID","SeqID","Abundance"))

# Load UC parquet file
cat("Loading UC membership table...\n")
UC <- fread(UC_FILE, sep="\t", header=FALSE)
UC <- unique(UC)
setnames(UC, c("SeqID","OTU"))

# Merge sequence table with cluster membership
SEQTAB <- merge(SEQTAB, UC, by="SeqID", all.x=TRUE)
SEQTAB <- SEQTAB[!is.na(OTU)]

cat("Unique sequences:", length(unique(SEQTAB$SeqID)), "\n")
cat("Clusters (OTUs):", length(unique(SEQTAB$OTU)), "\n")
cat("Samples:", length(unique(SEQTAB$SampleID)), "\n\n")

# Summarize abundance per OTU per sample
OTUTAB <- SEQTAB[, .(Abundance=sum(Abundance)), by=c("OTU","SampleID")]
OTUTAB[, Total := sum(Abundance), by="OTU"]

# UNCROSS tag-jump scoring
uncross_score <- function(x, N, n, f=0.01, tmin=0.1, p=1){
  z <- f * N / n
  sc <- 2 / (1 + exp(x/z)^p)
  data.table(Score=sc, TagJump=sc>=tmin)
}

cat("Calculating UNCROSS scores...\n")
OTUTAB <- cbind(OTUTAB, uncross_score(OTUTAB$Abundance, OTUTAB$Total, length(unique(OTUTAB$SampleID)), f=F_PARAM, p=P_PARAM))

cat("Tag-jump events:", sum(OTUTAB$TagJump), "\n\n")

# Export tag-jump scores
saveRDS(OTUTAB, file.path(OUTPUT_DIR,"TagJump_scores.rds"))

# Plot
PP <- ggplot(OTUTAB, aes(x=Total, y=Abundance, color=TagJump)) +
  geom_point() + scale_x_log10() + scale_y_log10() +
  scale_color_manual(values=c("#0C7C59","#D64933")) +
  labs(x="Total abundance of OTU (reads)", y="Abundance in sample (reads)")

pdf(file.path(OUTPUT_DIR,"TagJump_plot.pdf"), width=12, height=9.5, useDingbats=FALSE)
PP
dev.off()

# Calculate summary statistics
TJ <- data.table(
  Total_reads=sum(SEQTAB$Abundance),
  Number_of_TagJump_Events=sum(OTUTAB$TagJump),
  TagJump_reads=sum(OTUTAB[TagJump==TRUE]$Abundance)
)
TJ[, ReadPercent_removed := TagJump_reads / Total_reads * 100]
fwrite(TJ, file.path(OUTPUT_DIR,"TagJump_stats.txt"), sep="\t")

# Keep only non-tag-jump reads
RES <- merge(SEQTAB, OTUTAB[,.(OTU, SampleID, TagJump)], by=c("OTU","SampleID"))
RES <- RES[TagJump==FALSE, .(SampleID, SeqID, Abundance)]
setorder(RES, SampleID, -Abundance)

fwrite(RES, file.path(OUTPUT_DIR,"Seq_tab_TagJumpFiltered.txt"), sep="\t")

cat("All outputs written to", OUTPUT_DIR, "\n")

