import pandas as pd
import numpy as np
import progressbar as bar
import os
import scipy as sp
import yaml
import argparse
import matplotlib.pyplot as plt
import lifelines as lf
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr

parser = argparse.ArgumentParser()
parser.add_argument("--sig", dest="sig_file")
parser.add_argument("--outdir", dest="out_dir")

args = parser.parse_args()

with open(args.sig_file, 'r') as sigfile:
    signatures = yaml.safe_load(sigfile)

if "CANCER_TYPES" not in signatures:
    raise KeyError("Signature YAML file must contain a CANCER_TYPES key")

cancer_types = signatures.pop("CANCER_TYPES")
basedir = "/data/Robinson-SB/pancan-gc"

print("Reading clinical data")
clin = pd.read_excel(os.path.join(basedir, "PanCan_clin.xlsx"))

# Need to merge last_contact_days_to and death_days_to column
clin.last_contact_days_to = clin.last_contact_days_to.fillna(clin.death_days_to)

# Remove patients without any associated survival data
has_surv_data = ~(np.isnan(clin.death_days_to) & np.isnan(clin.last_contact_days_to))
clin = clin[has_surv_data]

print("Reading expression data")
exps = {ctype: pd.read_csv(os.path.join(basedir, "exps/%s_exp.tsv" % ctype), delim_whitespace=True) for ctype in
        cancer_types}

print("Reading purity data")
purity = pd.read_csv(os.path.join(basedir, "PanCan_purity.tsv"), delim_whitespace=True)
purity = purity[["array", "purity"]]
purity.array = [samp[:12] for samp in purity.array]

activity_pancan = pd.DataFrame(columns=exps.keys())

for signature in signatures:
    activity_pancan = activity_pancan.append(pd.Series(name=signature))

kmf = lf.KaplanMeierFitter()
mean_gene_exprs = {}
for ctype, exp in exps.items():
    exp.columns = [col[:12] for col in exp.columns]

    activity = pd.DataFrame(columns=exp.columns[1:])
    for signature in signatures:
        activity = activity.append(pd.Series(name=signature))

    # Only keep the gene symbols
    exp.gene_id = exp.gene_id.apply(lambda x: x.split("|")[0])
    exp = exp.loc[exp.gene_id != "?", :]

    # Get mean exprs of all genes in this cancer type
    print("[{0}] Getting mean gene expressions".format(ctype))
    mean_gene_exprs[ctype] = pd.DataFrame({"mean_exprs": exp.mean(axis=1), "delta_highctl": np.nan}).set_index(
        exp.gene_id)

    # Get GC activity score for each sample
    print("[{0}] Running signatures".format(ctype))
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

    corr_outdir = os.path.join(basedir, args.out_dir)
    if not os.path.exists(corr_outdir):
        os.mkdir(corr_outdir)

    corr_matrix.to_csv(os.path.join(corr_outdir, "%s_corr.tsv" % ctype), sep="\t")

    # Get mean expression across all samples in this cancer type
    for sig in activity.columns:
        activity_pancan[ctype][sig] = activity[sig].mean()

    # Build KM plot of this cancer type. Compare high-CTL patient survival against low-CTL
    clin_exprs = clin.merge(activity, left_on="bcr_patient_barcode", right_index=True).set_index("bcr_patient_barcode")
    clin_exprs = clin_exprs.merge(exp.set_index("gene_id").transpose(), left_index=True, right_index=True)
    high_ctl_activity_mask = (clin_exprs.ctl_activity > activity_pancan[ctype].ctl_activity)

    num_genes = len(mean_gene_exprs[ctype].index)
    sig_genes = []
    with bar.ProgressBar(max_value=num_genes, redirect_stdout=True) as gene_progress:
        i = 0
        for gene in mean_gene_exprs[ctype].index:
            gene = "TGFB1"
            i += 1
            gene_progress.update(i)

            # high_gene_mask = (clin_exprs[gene] > mean_gene_exprs[ctype]["mean_exprs"][gene])
            high_gene_mask = (sp.stats.zscore(clin_exprs[gene]) >= 1)
            high_ctl_high_gene = clin_exprs.loc[high_ctl_activity_mask & high_gene_mask]
            high_ctl_low_gene = clin_exprs.loc[high_ctl_activity_mask & ~]

            baseline_fit = kmf.fit(low_ctl_control.last_contact_days_to, event_observed=low_ctl_control.vital_status,
                                   label="Control")

            if len(high_ctl_high_gene) < 10:
                continue

            test_gene_fit = kmf.fit(high_ctl_high_gene.last_contact_days_to, event_observed=high_ctl_high_gene.vital_status,
                                    label="High CTL & %s" % gene)
            lrt = logrank_test(low_ctl_control.last_contact_days_to, high_ctl_high_gene.last_contact_days_to,
                               event_observed_A=low_ctl_control.vital_status,
                               event_observed_B=high_ctl_high_gene.vital_status)

            # Record genes that have a statistically significant difference vs baseline
            # TODO: FDR correct the p vals
            if lrt.p_value < 0.05:
                print("[%s: (lrt=%.2f, p=%.2f)]" % (gene, lrt.test_statistic, lrt.p_value))
                sig_genes.append(gene)



    # high_test_sig_mask = (sp.stats.zscore(clin_exprs.test_sig) >= 1)
    #
    # control_low_ctl = clin_exprs.loc[~high_ctl_activity_mask, :]
    # high_test_high_ctl = clin_exprs.loc[high_test_sig_mask & high_ctl_activity_mask, :]
    # low_test_high_ctl = clin_exprs.loc[~high_test_sig_mask & high_ctl_activity_mask, :]
    #
    # ax = plt.subplot(111)
    #
    # print(kmf.fit(high_test_high_ctl.last_contact_days_to, event_observed=high_test_high_ctl.vital_status, label="High test sig"))
    # ax = kmf.plot(ax=ax)
    # print(kmf.fit(control_low_ctl.last_contact_days_to, event_observed=control_low_ctl.vital_status, label="Control (low CTL)"))
    # ax = kmf.plot(ax=ax)
    # # plt.title("High GC response, low vs high CTL survival")
    # plt.show()
    #
    # ax = plt.subplot(111)
    # print(kmf.fit(low_test_high_ctl.last_contact_days_to, event_observed=low_test_high_ctl.vital_status, label="Low test sig"))
    # ax = kmf.plot(ax=ax)
    # print(kmf.fit(control_low_ctl.last_contact_days_to, event_observed=control_low_ctl.vital_status, label="Control (low CTL)"))
    # ax = kmf.plot(ax=ax)
    # # plt.title("GC")
    # plt.show()
    #
    # lrt = logrank_test(control_high_ctl.last_contact_days_to, low_gc_prod_high_ctl_act.last_contact_days_to,
    #                    event_observed_A=control_high_ctl.vital_status,
    #                    event_observed_B=low_gc_prod_high_ctl_act.vital_status)

activity_pancan = activity_pancan.transpose()
activity_pancan = activity_pancan.astype(float)
pancan_corrs = activity_pancan.corr(method="spearman")
pancan_corrs_pvals = pd.DataFrame(
    [[spearmanr(activity_pancan[rowsig], activity_pancan[colsig]).pvalue for rowsig in activity_pancan.columns] for
     colsig in activity_pancan.columns])
