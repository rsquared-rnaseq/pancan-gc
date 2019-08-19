# utils.py: contains simple python snippets I've used but don't warrant a whole separate file

import os
import glob
import pandas as pd


def merge_tsvs_by_cols(tsvs):
    the_tsv_df = pd.read_csv(tsvs[0], delim_whitespace=True)
    for tsv in tsvs[1:]:
        tsv_df = pd.read_csv(tsv, delim_whitespace=True)
        the_tsv_df = the_tsv_df.merge(tsv_df, left_on="Gene", right_on="Gene")

    return the_tsv_df


def simple_rbind(tsvs):
    tsvs_df = [pd.read_csv(tsv, delim_whitespace=True,
                           names=["barcode", "sample", "celltype", "treatment", "status", "uk1", "uk2"]) for tsv in tsvs]
    return pd.concat(tsvs_df, sort=False)


print("Getting expression data")
base_dir = "/data/Robinson-SB/pancan-gc/single_cell_datasets/GSE115978"
# cancer_outf = os.path.join(base_dir, "Mel_cancer_tpm_merged.tsv")
# noncancer_outf = os.path.join(base_dir, "Mel_noncancer_tpm_merged.tsv")
metadata_outf = os.path.join(base_dir, "meta_noncancer_cancer.tsv")

# cancer_tsvs = glob.glob(os.path.join(base_dir, "*_cancer_tpm.txt"))
# noncancer_tsvs = glob.glob(os.path.join(base_dir, "*_noncancer_tpm.txt"))
metadata = glob.glob(os.path.join(base_dir, "*_noncancer_cells.txt")) # Mel_53_noncancer_cells

# Eliminate empty metadata files
for meta_f in metadata:
    if os.stat(meta_f).st_size == 0:
        metadata.remove(meta_f)

print("Merging")
# cancer_merged = merge_tsvs_by_cols(cancer_tsvs)
# noncancer_merged = merge_tsvs_by_cols(noncancer_tsvs)
metadata_merged = simple_rbind(metadata)

# # Write out the merged files
print("Writing merged data")
# cancer_merged.to_csv(cancer_outf, sep='\t', index=False)
# noncancer_merged.to_csv(noncancer_outf, sep='\t', index=False)
# metadata_merged.to_csv(metadata_outf, sep='\t', index=False)
