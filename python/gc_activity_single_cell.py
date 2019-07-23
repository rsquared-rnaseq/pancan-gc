import pandas as pd
import os
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (8, 3)
base_dir = "/data/Robinson-SB/pancan-gc/single_cell_datasets/GSE115978"

# Start with glucocorticoid signatures
signatures = {
    "gc_production": ["STAR", "CYP11A1", "CYP11B1", "HSD11B1", "HSD11B2", "CYP21A2", "H6PD"],
    "gc_response": ["CD72", "DUSP1", "ETS1", "FKBP5", "FOS", "IL1RL1", "IRAK3", "IRF1", "MAP2K1", "MAP3K8", "NFKBIA",
                    "PER1", "PIK3CA", "PIK3R1", "SGK1", "THBD", "TLR2", "TNFAIP3", "TSC22D3"],
    "cct3": ["CCT3"]
}
signature_out = pd.DataFrame(columns=["cancer", "noncancer"])

print("Reading in data... ", end="")
metadata = pd.read_csv(os.path.join(base_dir, "Mel_metadata_merged.tsv"), delim_whitespace=True)
cancer_tpm = pd.read_pickle(os.path.join(base_dir, "Mel_cancer_tpm_merged.pkl"))
noncancer_tpm = pd.read_pickle(os.path.join(base_dir, "Mel_noncancer_tpm_merged.pkl"))
print("done")

for sig_name in signatures.keys():
    signature_out = signature_out.append(pd.Series(name=sig_name))

for gc_prod_gene in signatures["cct3"]:
    plt.plot()
    # Different cell types in the noncancer sample
    types = ['B.cell', 'CAF', 'Endo.', 'Macrophage', 'NK', 'T.CD4', 'T.CD8', 'T.cell']

    data_plot = []  # noncancer_tpm[gc_prod_gene + "_NC"]
    plot_labels = []

    cancer_tpm_withzeros = cancer_tpm[gc_prod_gene + "_C"]
    cancer_tpm_nozeros = cancer_tpm_withzeros.loc[cancer_tpm_withzeros != 0]
    percent_gene_in_cancer = len(cancer_tpm_nozeros) / len(cancer_tpm_withzeros) * 100
    cancer_plotlabel = "Cancer\n{0}/{1}\n({2:.1f}%)".format(len(cancer_tpm_nozeros), len(cancer_tpm_withzeros),
                                                            percent_gene_in_cancer)
    plot_labels.append(cancer_plotlabel)
    data_plot.append(cancer_tpm_withzeros)
    for celltype in types:
        cells_with_type = metadata.loc[metadata.celltype == celltype].barcode
        tpms_type = noncancer_tpm.loc[cells_with_type, :]
        tpms_type_with_gene_withzeros = tpms_type[gc_prod_gene + "_NC"]
        tpms_type_with_gene_nozeros = tpms_type_with_gene_withzeros.loc[tpms_type_with_gene_withzeros != 0]

        data_plot.append(tpms_type_with_gene_withzeros)
        percent_gene_in_ctype = len(tpms_type_with_gene_nozeros) / len(tpms_type_with_gene_withzeros) * 100
        plot_labels.append(
            "{0}\n{1}/{2}\n({3:.1f}%)".format(celltype, len(tpms_type_with_gene_nozeros),
                                              len(tpms_type_with_gene_withzeros), percent_gene_in_ctype))

    plt.boxplot(data_plot, labels=plot_labels)
    plt.title(gc_prod_gene)
    plt.show()
