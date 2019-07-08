import pandas as pd
import numpy as np
import progressbar as bar
import os
import scipy as sp
import yaml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sig", dest="sig_file")
parser.add_argument("--outdir", dest="out_dir")

args = parser.parse_args()

with open(args.sig_file, 'r') as sigfile:
    signatures = yaml.safe_load(sigfile)

cancer_types = ["SKCM"]
# cancer_types = ["BRCA", "GBM", "OV", "LUAD", "UCEC", "KIRC", "HNSC", "LGG", "THCA", "LUSC", "PRAD", "SKCM", "COAD",
#                 "STAD", "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML", "ESCA", "PAAD", "PCPG", "READ", "TGCT", "THYM",
#                 "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS", "CHOL"]
basedir = "/data/Robinson-SB/pancan-gc"

# print("Reading clinical data")
# clin = pd.read_excel(os.path.join(basedir, "PanCan_clin.xlsx"))

print("Reading expression data")
exps = {ctype: pd.read_csv(os.path.join(basedir, "exps/%s_exp.tsv" % ctype), delim_whitespace=True) for ctype in cancer_types}
# TODO: Support PanCan
# exp = pd.concat(exps, ignore_index=True)

print("Reading purity data")
purity = pd.read_csv(os.path.join(basedir, "PanCan_purity.tsv"), delim_whitespace=True)
purity = purity[["array", "purity"]]
purity.array = [samp[:12] for samp in purity.array]

for ctype, exp in exps.items():
    exp.columns = [col[:12] for col in exp.columns]

    # Set up DataFrame to store our results
    activity = pd.DataFrame(columns=exp.columns[1:])
    for signature in signatures:
        activity = activity.append(pd.Series(name=signature))

    # Only keep the gene symbols
    exp.gene_id = exp.gene_id.apply(lambda x: x.split("|")[0])

    # Get GC activity score for each sample
    print("Running signatures for", ctype)
    with bar.ProgressBar(max_value=len(exp.columns[1:])) as progress:
        i = 0
        for sample in exp.columns[1:]:
            # For each sample, compute all the signatures by mean expression value
            for sig_name, signature in signatures.items():
                activity.loc[sig_name][sample] = exp.loc[exp.gene_id.isin(signature), sample].mean()

            i += 1
            progress.update(i)

    # Convert the columns to floats. This allows us to calculate correlations as well as makes computations much faster
    activity = activity.astype(float)

    # Transpose the DataFrame. This allows us to use the corr() method to calculate correlations between columns,
    # and makes the DataFrame look nicer in SciView
    # Merge purity and sample data
    activity = activity.transpose().merge(purity, left_index=True, right_on="array", how="outer")
    activity = activity[activity.array.isin(exp.columns)]
    activity = activity[pd.notnull(activity.purity)]
    activity = activity.set_index("array")

    corr_matrix = activity.corr(method="spearman")
    corr_matrix.to_csv(os.path.join(basedir, "%s/%s_corr.tsv" % (args.out_dir, ctype)), sep="\t")
