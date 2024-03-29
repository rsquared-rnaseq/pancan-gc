import pandas as pd
import numpy as np
import progressbar as bar
import os
import yaml
import lifelines as lf
from lifelines.statistics import logrank_test
import scipy.stats as ss
import matplotlib.pyplot as plt

sig_file = "../signatures/tumor_purity_analysis.yml"
out_dir = "corrs_liver"


def get_pvals(r, n):
    # uses t-test to get p-values
    t = r * np.sqrt((n - 2) / (1 - r * r))
    return ss.t.cdf(t, n - 2)


with open(sig_file, 'r') as sigfile:
    signatures = yaml.safe_load(sigfile)

if "CANCER_TYPES" not in signatures:
    raise KeyError("Signature YAML file must contain a CANCER_TYPES key")

cancer_types = signatures.pop("CANCER_TYPES")
basedir = "/data/Robinson-SB/pancan-gc"

print("Reading clinical data")
clin = pd.read_pickle(os.path.join(basedir, "PanCan_clin.pkl"))

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
    pval_matrix = get_pvals(corr_matrix, len(corr_matrix))

    corr_outdir = os.path.join(basedir, out_dir)
    if not os.path.exists(corr_outdir):
        os.mkdir(corr_outdir)

    corr_matrix.to_csv(os.path.join(corr_outdir, "%s_corr.tsv" % ctype), sep="\t")

    # Get mean expression across all samples in this cancer type
    for sig in activity.columns:
        activity_pancan[ctype][sig] = activity[sig].mean()

    # Build KM plot of this cancer type. Compare high-CTL patient survival against low-CTL
    clin_exprs = clin.merge(activity, left_on="bcr_patient_barcode", right_index=True).set_index("bcr_patient_barcode")
    clin_exprs = clin_exprs.merge(exp.set_index("gene_id").transpose(), left_index=True, right_index=True)

    # high_cct3_mask = (clin_exprs.cct3 > activity_pancan[ctype].cct3)
    top_quartile_cct3 = (clin_exprs.cct3 > np.percentile(clin_exprs.cct3, q=75))
    bottom_quartile_cct3 = (clin_exprs.cct3 < np.percentile(clin_exprs.cct3, q=25))

    high_cct3 = clin_exprs.loc[top_quartile_cct3]
    low_cct3 = clin_exprs.loc[bottom_quartile_cct3]

    ax = plt.subplot(111)
    km5 = kmf.fit(high_cct3.last_contact_days_to, event_observed=high_cct3.vital_status, label="Top quartile CCT3")
    ax = kmf.plot(ax=ax)
    km6 = kmf.fit(low_cct3.last_contact_days_to, event_observed=low_cct3.vital_status, label="Bottom quartile CCT3")
    ax = kmf.plot(ax=ax)

    plt.title("Top vs bottom quartile CCT3 expression (%s)" % ctype)
    plt.show()

    # Calculate statistical significance
    lrt = logrank_test(high_cct3.last_contact_days_to, low_cct3.last_contact_days_to,
                       event_observed_A=high_cct3.vital_status, event_observed_B=low_cct3.vital_status)
    #
    # high_ctl_activity_mask = (clin_exprs.ctl_activity > activity_pancan[ctype].ctl_activity)
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

activity_pancan = activity_pancan.transpose()
activity_pancan = activity_pancan.astype(float)
pancan_corrs = activity_pancan.corr(method="spearman")
