#!/usr/bin/env python
import sys
import argparse
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('Agg')
from pathlib import Path
import matplotlib.pyplot as plt



def get_ase(ase_fname, n=None):
    """parse the ASE TSV for the top n most promising str-tissue pairs"""
    ase = pd.read_csv(
        ase_fname, sep="\t", header=0, index_col=['str', 'tissue'], usecols=[
            'str', 'tissue', 'chrom', 'str_pos', 'str_a1', 'str_a2', 'distance',
            'variant_annotation', 'gene', 'ase', 'allelic_ratio'
        ]
    )
    str_tissue_pairs = ase.index.unique()
    if n is None:
        n = len(str_tissue_pairs)
    # get the top n most promising str-tissue pairs
    ase = ase.loc[str_tissue_pairs[:n],:]
    ase['repeat_ratio'] = ase['str_a1'].str.len()/ase['str_a2'].str.len()
    return ase


def plot(ase, out):
    """make plots for each str-tissue pair in the ase table"""
    for STR in ase.index.unique():
        plt.figure()
        repeat = ase['repeat_ratio'][STR].to_numpy()
        allelic = ase['allelic_ratio'][STR].to_numpy()
        # perform a simple linear regression
        p = np.polyfit(repeat, allelic, 1)
        r_squared = 1 - (
            sum(
                (allelic - (p[0] * repeat + p[1]))**2
            ) / (
                (len(allelic) - 1) * np.var(allelic, ddof=1)
            )
        )
        p = np.poly1d(p)
        # plot the points and the line
        plt.scatter(repeat, allelic, color='r', label="_nolegend_")
        plt.xlabel("Ratio of STR Lengths")
        plt.ylabel("Allelic Expression Ratio")
        plt.plot(
            repeat,
            p(repeat),
            label=str(p)+"\nr-squared: "+str(round(r_squared, 2))
        )
        # TODO: include information about distance, tissue, gene?, variant_annotation?
        plt.legend(frameon=False, loc='lower right')
        plt.tight_layout()
        plt.savefig(str(out)+"/"+"-".join(STR)+".png", bbox_inches='tight')
        plt.clf()



def get_args():
    """parse the arguments to this script using argparse"""
    parser = argparse.ArgumentParser(description=
        "Plot linear associations between ASE and STR repeat number."
    )
    parser.add_argument(
        '-n', '--num-plots', default=None, type=int,
        help='Only create plots for the first n STRs'
    )    
    parser.add_argument(
        'out', type=Path, help='Path to a directory in which to place each plot'
        # TODO: describe the structure of this tsv file better
    )
    parser.add_argument(
        'ase', default=sys.stdin, nargs='?',
        help="""
            A TSV containing ASE for each STR. The columns of this TSV must have names:
            str, tissue, chrom, str_pos, str_a1, str_a2, distance, variant_annotation,
            gene, ase
            If this argument is not provided, it will be read from stdin
        """
    )
    args = parser.parse_args()
    return args


def main(args):
    # create the output directory if it doesn't exist yet
    args.out.mkdir(exist_ok=True)
    plot(get_ase(args.ase, args.num_plots), args.out)


if __name__ == "__main__":
    main(get_args())