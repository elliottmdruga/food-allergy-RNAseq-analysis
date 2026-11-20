"""
Analysis pipeline following RNAseq quantification with nf-core
Should work with any usual output from the star-salmon pipeline
"""

# Imports
import sys
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# want matplotlib for volcano plots
# want seaborn clustermap for heatmaps

"""----------INPUT PROCESSING----------"""

# Arguments that will print the overview and usage message, and those messages
help_strings = ["h", "help", "-h", "-help", "--h", "--help"]
overview_message = ("Script for producing downstream analyses of RNAseq data processed by nf-core"
                    " RNAseq pipeline with star-salmon. Generates a variety of visualizations.")
usage_message = "USAGE: python analysis.py [counts tsv file] [metadata tsv file]"

# Check if the first argument provided is a help string, print usage if so
if len(sys.argv) > 1 and sys.argv[1] in help_strings:
    print(overview_message)
    print(usage_message)
    sys.exit()

# Check if enough arguments are provided, and print a message if not
if len(sys.argv) < 3:
    print("Insufficient arguments.")
    print(usage_message)
    sys.exit()

# Read filenames for counts and metadata files from arguments
counts_file = sys.argv[1]
metadata_file = sys.argv[2]

#TODO - make these optional arguments to the program

# Low-expression cutoff for filtering out genes
count_cutoff = 100
# Whether to save results and the file path
save_results = False
output_path = "./"

# Read counts file into a dataframe
try:
    counts = pd.read_csv(counts_file, sep="\t")
# Print a message if file is not found
except FileNotFoundError:
    print(f"{counts_file} not found. Please provide a valid counts tsv file.")
    print(usage_message)
    sys.exit()

# Read metadata file into a dataframe
try:
    meta = pd.read_csv(metadata_file, index_col=0)
# Print a message if file is not found
except FileNotFoundError:
    print(f"{metadata_file} not found. Please provide a valid counts csv file.")
    print(usage_message)
    sys.exit()

"""----------FILTERING AND PREPROCESSING----------"""

# Get rid of gene_id column - not needed
counts = counts.drop(columns=["gene_id"])

# Set gene_id and gene_name as index columns
counts = counts.set_index(['gene_name'])

# Convert raw counts to integers for analysis in Deseq2
counts = counts.round().astype(int)

# Filter out low-count genes according to cutoff
counts = counts.loc[counts.sum(axis=1) >= count_cutoff, :]

# Transpose counts matrix for use in Deseq2
counts = counts.T
# Rename index column
counts.index.name = 'sample'

# Ensure index column is the same for both counts and metadata
meta = meta.loc[counts.index]

sns.clustermap(counts)
plt.show()

"""----------ANALYSIS----------"""
"""
# Create deseq dataset from counts and metadata
# Using group study design - no other significant varaibles to consider
# Maybe make this more variable based on input?
dds = DeseqDataSet(counts=counts, metadata=meta, design="~ group")
dds.deseq2()
stats = DeseqStats(dds, contrast=["group", "FA_NA", "FA_PN"])
stats.summary()
res = stats.results_df
"""
"""----------VISUALIZATION----------"""

