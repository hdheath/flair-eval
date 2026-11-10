#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import sys

file_path = "/private/groups/brookslab/hdheath/projects/flair-eval/results/ted/partition_test/A549_cDNA/A549_cDNA_default_all_collapse_with_gtf_default_collapse_ted.tsv"
df = pd.read_csv(file_path, sep='\t')

metrics = ['5prime_precision', '5prime_recall', '3prime_precision', '3prime_recall']
values = [df[metric].iloc[0] for metric in metrics]

# Create bar plot
plt.figure(figsize=(10, 6))
plt.bar(metrics, values, color=['blue', 'green', 'orange', 'red'])
plt.xticks(rotation=45)
plt.ylabel('Values')
plt.title('TED Results Summary')
plt.tight_layout()
plt.savefig('ted_results_plot.png')
plt.show()
