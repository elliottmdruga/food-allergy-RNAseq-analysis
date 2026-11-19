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

# Low-expression cutoff for filtering out genes
count_cutoff = 10000 #TODO - add as an optional command line parameter

# Set gene_id and gene_name as index columns
counts = counts.set_index(['gene_id', 'gene_name'])

# Convert raw counts to integers for analysis in Deseq2
counts = counts.round().astype(int)

# Filter out low-count genes according to cutoff
counts = counts.loc[counts.sum(axis=1) >= count_cutoff, :]

print(counts.head())

"""----------ANALYSIS----------"""
"""
dds = DeseqDataSet(counts=counts, metadata=meta, design="~ sex + group")
dds.deseq2()
stats = DeseqStats(dds, contrast=["group", "FA", "NA"])
stats.summary()
res = stats.results_df
"""
"""----------VISUALIZATION----------"""