import pandas as pd
import numpy as np
import progressbar as bar
import os
from scipy.stats import spearmanr

# TODO (gc_activity sig): Matt's email mentions CYP21A1, but this gene only exists in mice. In humans the same gene
# TODO: (21-hydroxylase) is encoded by the CYP21A2 gene. Which one should I use?
signatures = {
    # "gc_activity": ["STAR", "CYP11A1", "CYP11B1", "HSD11B1", "HSD11B2", "CYP21A2"],
    "gc_activity": ["STAR", "CYP11A1", "CYP11B1", "HSD11B1", "HSD11B2", "CYP21A1", "HSD11B1"],
    "ctl_activity": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1"]
}
cancer_type = "SKCM"
basedir = "/data/Robinson-SB/pancan-gc"

print("Reading clinical data")
clin = pd.read_excel(os.path.join(basedir, "PanCan_clin.xlsx"))

print("Reading expression data")
exp = pd.read_csv(os.path.join(basedir, "exps/%s_exp.tsv" % cancer_type), delim_whitespace=True)

# Set up DataFrame to store our results
activity = pd.DataFrame(columns=exp.columns[1:])

activity = activity.append(pd.Series(name="gc_activity"))
activity = activity.append(pd.Series(name="ctl_activity"))

# Only keep the gene symbols
exp['gene_id'] = exp['gene_id'].apply(lambda x: x.split("|")[0])
# exp.loc[exp['gene_id'].str.contains("A1BG")]

# Get GC activity score for each sample
print("Running signatures")
with bar.ProgressBar(max_value=len(exp.columns[1:])) as progress:
    i = 0
    for sample in exp.columns[1:]:
        # For each sample, compute all the signatures by mean expression value
        for sig_name, signature in signatures.items():
            activity.loc[sig_name][sample] = exp.loc[exp["gene_id"].isin(signature), sample].mean()
            # samp_sig = 0
            # for sig_gene in signature:
            #     sig_gene_exp = exp.loc[exp["gene_id"] == sig_gene][sample]
            #     if len(sig_gene_exp) == 1:
            #         samp_sig += sig_gene_exp.item()
            #     else:
            #         print("disparity gene=%s, len=%d" % (sig_gene, len(sig_gene_exp)))
            #         print(sig_gene_exp)
            # samp_sig /= len(signature)
            # activity.loc[sig_name][sample] = samp_sig

        i += 1
        progress.update(i)

# Transpose the DataFrame. This allows us to use the corr() method to calculate correlations between columns,
# and makes the DataFrame look nicer in SciView
activity = activity.transpose()

# Convert the columns to floats. This allows us to calculate correlations as well as makes computations much faster
activity = activity.astype(float)



