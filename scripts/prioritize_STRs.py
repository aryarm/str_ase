#!/usr/bin/env python
import sys
import pandas as pd



# how should we weight the values that contribute to the importance metric?
WEIGHTS = {
    'num_samples': 6,
    'allele_variance': 5,
    'distance_mean': 1,
    'ase_mean': 4,
    'ase_variance': 7
}

infile = sys.argv[1] if len(sys.argv) >= 2 else sys.stdin
outfile = sys.argv[2] if len(sys.argv) >= 3 else sys.stdout
counts = pd.read_csv(infile, sep="\t", header=0, index_col=['str', 'tissue', 'subject'])

def get_ase(counts):
    """defined as abs(0.5-(a1_count/total_count))"""
    return (0.5-counts['a1_count']/counts['total_count']).abs()

def normalize(vals):
    """normalize a series of values using min/max normalization"""
    return (vals-vals.min())/(vals.max()-vals.min())

# normalize the weights
weight_sum = sum(WEIGHTS.values())
WEIGHTS = {metric: WEIGHTS[metric]/weight_sum for metric in WEIGHTS}

# calculate all of our metrics
ase = get_ase(counts[['a1_count', 'total_count']])
metrics = {
    'num_samples': counts.groupby(['str', 'tissue']).size(),
    'allele_variance': counts[['str_a1', 'str_a2']].stack().groupby(
        ['str', 'tissue']
    ).apply(lambda x: x.drop_duplicates().str.len().std()),
    'distance_mean': counts['distance'].groupby(['str', 'tissue']).mean(),
    'ase_mean': ase.groupby(['str', 'tissue']).mean(),
    'ase_variance': ase.groupby(['str', 'tissue']).std().fillna(0)
}
# normalize each metric
metrics = {metric: normalize(metrics[metric]) for metric in metrics}

assert set(metrics.keys()) == set(WEIGHTS.keys())
# the importance is the weighted sum of each metric
counts['importance'] = sum([WEIGHTS[metric]*metrics[metric] for metric in metrics])

counts.sort_values('importance', ascending=False).to_csv(outfile, sep="\t")
