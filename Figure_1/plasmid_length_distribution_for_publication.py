# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 09:39:56 2025

@author: hayat
"""

import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
import seaborn as sns
import pandas as pd

plt.rcParams["svg.fonttype"] = "none"  
plt.rcParams["font.family"] = "Arial"  
plt.rcParams["font.size"] = 10

def analyze_fasta(fasta_file):
    # Read sequences and calculate lengths
    lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    
    # Calculate statistics
    min_length = min(lengths)
    max_length = max(lengths)
    median_length = np.median(lengths)
    
    # Print statistics
    print(f"Minimum length: {min_length}")
    print(f"Maximum length: {max_length}")
    print(f"Median length: {median_length}")
    
    # Plot histogram
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=50, edgecolor='black')
    plt.title("Distribution of Sequence Lengths for Top Widespread and Abundant Plasmids")
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    
    # Add vertical lines for min, max, and median
    plt.axvline(min_length, color='r', linestyle='dashed', linewidth=2, label=f'Min: {min_length}')
    plt.axvline(max_length, color='g', linestyle='dashed', linewidth=2, label=f'Max: {max_length}')
    plt.axvline(median_length, color='b', linestyle='dashed', linewidth=2, label=f'Median: {median_length}')
    
    plt.legend()
    #plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/sequence_length_distribution_all_sponge_plas.png", dpi=300)
    plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/sequence_length_distribution_all_sponge_plas.svg")
    plt.close()

#analyze_fasta('C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids.fasta')




def analyze_fasta_kde(fasta_file):
    lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

    min_length = min(lengths)
    max_length = max(lengths)
    median_length = np.median(lengths)

    print(f"Minimum length: {min_length}")
    print(f"Maximum length: {max_length}")
    print(f"Median length: {median_length}")

    plt.figure(figsize=(10, 6))
    sns.kdeplot(lengths, fill=True, color="skyblue", bw_adjust=0.5)
    #plt.title("Density of Sequence Lengths (Plasmids)")
    plt.xlabel("Sequence Length")
    plt.ylabel("Density")

    plt.axvline(min_length, color='r', linestyle='dashed', linewidth=2, label=f'Min: {min_length}')
    plt.axvline(max_length, color='g', linestyle='dashed', linewidth=2, label=f'Max: {max_length}')
    plt.axvline(median_length, color='b', linestyle='dashed', linewidth=2, label=f'Median: {median_length}')

    plt.legend()
    plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/sequence_length_density.svg")
    plt.show()
    plt.close()

#analyze_fasta_kde('C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids.fasta')

def analyze_fasta_kde_log(fasta_file):
    lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

    min_length = min(lengths)
    max_length = max(lengths)
    median_length = np.median(lengths)

    print(f"Minimum length: {min_length}")
    print(f"Maximum length: {max_length}")
    print(f"Median length: {median_length}")
    df = pd.DataFrame({"length": lengths})
    df.to_csv("C:/Users/hayat/Downloads/R_files/data/widespread_plasmid_lengths.csv",
          index=False)

    plt.figure(figsize=(10, 6))
    sns.histplot(lengths, bins=50, kde=True, color="skyblue", log_scale=True)
    #plt.title("Distribution of Sequence Lengths (Histogram + KDE, Log Scale)")
    plt.xlabel("Sequence Length (log scale)")
    plt.ylabel("Count")


    plt.axvline(min_length, color='r', linestyle='dashed', linewidth=2, label=f'Min: {min_length}')
    plt.axvline(max_length, color='g', linestyle='dashed', linewidth=2, label=f'Max: {max_length}')
    plt.axvline(median_length, color='b', linestyle='dashed', linewidth=2, label=f'Median: {median_length}')

    plt.legend()
    
    #plt.savefig("C:/Users/hayat/Downloads/R_files/graphs/allplas_sequence_length_density_log.svg")
    plt.show()
    plt.close()

#analyze_fasta_kde_log('C:/Users/hayat/Downloads/R_files/data/all_sponge_plasmids.fasta')

analyze_fasta_kde_log('C:/Users/hayat/Downloads/R_files/data/top_abundant_and_widespread_final_plasmid.fasta')

