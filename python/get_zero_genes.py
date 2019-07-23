import pandas as pd
import numpy as np
import progressbar
import os

base_dir = "/data/Robinson-SB/pancan-gc/single_cell_datasets/GSE115978"

print("Reading in data... ", end="")
metadata = pd.read_csv(os.path.join(base_dir, "Mel_metadata_merged.tsv"), delim_whitespace=True)
cancer_tpm = pd.read_pickle(os.path.join(base_dir, "Mel_cancer_tpm_merged.pkl"))
noncancer_tpm = pd.read_pickle(os.path.join(base_dir, "Mel_noncancer_tpm_merged.pkl"))
print("done")

cancer_tpm.columns = [col[:-2] for col in cancer_tpm.columns]
noncancer_tpm.columns = [col[:-3] for col in noncancer_tpm.columns]

print("Getting cancer-only genes")
# Get genes that are only expressed in melanoma and no other cell types
mel_only = []
with progressbar.ProgressBar(max_value=len(cancer_tpm.columns)) as bar:
    i = 0
    for gene in cancer_tpm.columns:
        if cancer_tpm[gene].sum() > 0 and noncancer_tpm[gene].sum() == 0:
            mel_only.append(gene)

        i += 1
        bar.update(i)

print("Getting T-cell only genes")
tcell_samples = metadata.loc[metadata.celltype.str.startswith("T")].barcode
tcell_tpm = noncancer_tpm.loc[tcell_samples, :]

# Get genes that are only expressed in non-cancer cells and not cancer
tcell_only = []
with progressbar.ProgressBar(max_value=len(tcell_tpm.columns)) as bar:
    i = 0
    for gene in tcell_tpm.columns:
        if tcell_tpm[gene].sum() > 0 and cancer_tpm[gene].sum() == 0:
            tcell_only.append(gene)

        i += 1
        bar.update(i)