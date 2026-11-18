"""Analysis pipeline following RNAseq quantification with nf-core"""
"""Should work with any usual output from the star-salmon pipeline"""

import sys
import numpy as np
import pandas as pd

counts_file = sys.argv[1]
metadata_file = sys.argv[2]
help_strings = ["h", "help", "-h", "-help", "--h", "--help"]

counts = pd.read_csv(counts_file, sep="\t", index_col=0)
gene_names = counts.pop("gene_name")

print(counts.head())