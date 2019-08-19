import numpy as np
import pandas as pd
import os

base_dir = "/data/Robinson-SB/pancan-gc/checkpoint_datasets"
dat = pd.read_csv(os.path.join(base_dir, "van-allen_mRNA.tsv"), delim_whitespace=True)