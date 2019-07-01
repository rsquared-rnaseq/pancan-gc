configfile: "config.yml"

CANCERS = ["BRCA", "GBM", "OV", "LUAD", "UCEC", "KIRC", "HNSC", "LGG", "THCA", "LUSC", "PRAD", "SKCM", "COAD",
           "STAD", "BLCA", "LIHC", "CESC", "KIRP", "SARC", "LAML", "ESCA", "PAAD", "PCPG", "READ", "TGCT", "THYM",
           "KICH", "ACC", "MESO", "UVM", "DLBC", "UCS", "CHOL"]

EXPRESSIONF = config["pancan_dir"] + "/PanCan_exp.tsv"
CLINF = config["pancan_dir"] + "/PanCan_clin.xlsx"


rule all:
    input: expand(config["pancan_dir"] + "/exps/{cancer}_exp.tsv", cancer=CANCERS)


rule dl:
    input: EXPRESSIONF, CLINF

rule dl_pancan_exp:
    output: exp_out=EXPRESSIONF
    shell: "wget {config[pancan_exp_dl]} -O {output.exp_out}"

rule dl_pancan_clin:
    output: clin_out=CLINF
    shell: "wget {config[pancan_clin_dl]} -O {output.clin_out}"

rule split_exp:
    input: exp=EXPRESSIONF, clin=CLINF
    output:
        expand(config["pancan_dir"] + "/exps/{cancer}_exp.tsv", cancer=CANCERS)
    params:
        out_dir=config["pancan_dir"] + "/exps"
    script:
        "R/get_cancertype_exp.R"
