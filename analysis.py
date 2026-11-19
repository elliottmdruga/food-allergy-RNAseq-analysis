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


# Arguments that will print the usage message, and the message itself
help_strings = ["h", "help", "-h", "-help", "--h", "--help"]
usage_message = "USAGE: python analysis.py [counts tsv file] [metadata tsv file]"

# Check if the first argument provided is a help string, print usage if so
if len(sys.argv) > 1 and sys.argv[1] in help_strings:
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
    counts = pd.read_csv(counts_file, sep="\t", index_col=0)
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

print(counts.head())