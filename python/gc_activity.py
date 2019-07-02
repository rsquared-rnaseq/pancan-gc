import pandas as pd
import numpy as np
import progressbar as bar
import os
import scipy as sp

# TODO (gc_activity sig): Matt's email mentions CYP21A1, but this gene only exists in mice. In humans the same gene
# TODO: (21-hydroxylase) is encoded by the CYP21A2 gene. Which one should I use?
signatures = {
    "gc_production": ["STAR", "CYP11A1", "CYP11B1", "HSD11B1", "HSD11B2", "CYP21A1", "HSD11B1"],
    "gc_response": ["CD72", "DUSP1", "ETS1", "FKBP5", "FOS", "IL1RL1", "IRAK3", "IRF1", "MAP2K1", "MAP3K8", "NFKBIA",
                    "PER1", "PIK3CA", "PIK3R1", "SGK1", "THBD", "TLR2", "TNFAIP3", "TSC22D3"],
    "ctl_activity": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1"],
    "ctl_exhaustion": ["TNFSF6", "PBX3", "GP49B", "CD244", "CCL3", "EOMES", "CASP3", "PLSCR1", "KDT1", "CTLA4", "PDCD1",
                       "IER5", "RGS16", "A430109M19RIK", "TNFRSF9", "PENK1", "EOMES", "COCH", "PTPN13", "TCRG-V4",
                       "NR4A2", "CD160", "PTGER4", "CCL4", "WBP5", "GPR56", "1110067D22RIK", "ENTPD1", "SH2D2A",
                       "SEPT4", "ISG20", "TRIM47", "SERPINA3G", "CASP4", "9130009C22RIK", "C79248", "LAG3", "NR4A2",
                       "NFTAC1", "CAR2", "C330007P06RIK", "GPD2", "2700084L22RIK", "RNF11", "CAPZB", "TUBB2", "BUB1",
                       "JAK3", "9130410M22RIK", "CD9", "TCRG-V4", "1810054D07RIK", "RCN", "2010100O12RIK", "SYBL1",
                       "ETF1", "CPA3", "CD7", "ART3", "1810035L17RIK", "ATF1", "PRKWNK1", "MTV43", "CIT", "CCRL2",
                       "ADFP", "D8ERTD531E", "TCEA2", "MYH4", "TNFRSF1A", "SPP1", "S100A13", "PON2", "AI181996", "G1P2",
                       "TANK", "SHKBP1", "2510004L0RIK", "D15ERTD781E", "ICSBP1", "BC024955", "GDF3", "ITGAV",
                       "1110006I15RIK", "CPSF2", "KLK6", "CPT2", "LMAN2", "TOR3A", "CRYGB", "GPR65", "MKI67",
                       "TCRB-V13", "NPTXR", "SNRPB2", "NDFIP1", "PTGER2", "ZFP91", "SPOCK2", "5730469M10RIK", "CXCL10",
                       "GCIN", "TRIM25", "WBSCR5", "MOX2", "DOCK7", "PAWR", "CHL1"]
}
cancer_types = ["SKCM"]
basedir = "/data/Robinson-SB/pancan-gc"

# print("Reading clinical data")
# clin = pd.read_excel(os.path.join(basedir, "PanCan_clin.xlsx"))

print("Reading expression data")
exps = [pd.read_csv(os.path.join(basedir, "exps/%s_exp.tsv" % ctype), delim_whitespace=True) for ctype in cancer_types]
exp = pd.concat(exps, ignore_index=True)

# Set up DataFrame to store our results
activity = pd.DataFrame(columns=exp.columns[1:])

for signature in signatures:
    activity = activity.append(pd.Series(name=signature))

# Only keep the gene symbols
exp["gene_id"] = exp["gene_id"].apply(lambda x: x.split("|")[0])
# exp.loc[exp['gene_id'].str.contains("A1BG")]

# Get GC activity score for each sample
print("Running signatures")
with bar.ProgressBar(max_value=len(exp.columns[1:])) as progress:
    i = 0
    for sample in exp.columns[1:]:
        # For each sample, compute all the signatures by mean expression value
        for sig_name, signature in signatures.items():
            activity.loc[sig_name][sample] = exp.loc[exp["gene_id"].isin(signature), sample].mean()

        i += 1
        progress.update(i)

# Transpose the DataFrame. This allows us to use the corr() method to calculate correlations between columns,
# and makes the DataFrame look nicer in SciView
activity = activity.transpose()

# Convert the columns to floats. This allows us to calculate correlations as well as makes computations much faster
activity = activity.astype(float)

corr_matrix = activity.corr(method="spearman")


