library(data.table)
library(openxlsx)
library(dplyr)

# expressionFn <- snakemake@input[["exp"]]
# clinFn <- snakemake@input[["clin"]]
# outDir <- snakemake@output[["out_dir"]]
expressionFn <- "/data/Robinson-SB/pancan-gc/PanCan_exp.tsv"
clinFn <- "/data/Robinson-SB/pancan-gc/PanCan_clin.xlsx"
outDir <- "/data/Robinson-SB/pancan-gc/exps"

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

# select_cancers = unique(typePatients$bcr_patient_barcode)
# colnames(exp)
# 
# select_cancers = c(select_cancers, "gene_id")
# 
# exp_subset = exp[,.SD,.SDcols=c(select_cancers)]
# 
# select_cancers[1:10]
# colnames(exp)[1:10]
# 
# bob = intersect(colnames(exp), select_cancers)
# exp_subset = exp[,.SD,.SDcols=bob]
# 
# class(typePatients)
# whale = as.data.table(typePatients)
# type_subset = whale[bcr_patient_barcode %in% colnames(exp),]
# 
# type_subset = unique(type_subset)
# length(unique(type_subset$bcr_patient_barcode))
# 
# turtle = unique(select_cancers)
