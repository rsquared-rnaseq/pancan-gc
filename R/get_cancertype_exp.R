library(data.table)
library(openxlsx)
library(dplyr)

expressionFn <- snakemake@input[["exp"]]
clinFn <- snakemake@input[["clin"]]
outDir <- snakemake@params[["out_dir"]]

message("Reading in expression data")
exp <- fread(expressionFn)
clin <- read.xlsx(clinFn, sheet = 1, colNames = TRUE)
clin <- as.data.table(clin)

cancerTypes <- unique(clin$type)

for (ctype in cancerTypes) {
  typePatients <- subset(clin, type == ctype)
  expTypeSelector <- paste0("^(gene_id|", paste(typePatients$bcr_patient_barcode, collapse = "|"), ")")
  expType <- select(exp, matches(expTypeSelector))
  message(paste0("Writing TSV for cancer type ", ctype))
  fwrite(expType, paste0(outDir, "/", ctype, "_exp.tsv"), sep = "\t")
}
