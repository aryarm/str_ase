#!/usr/bin/env python
import sys
import csv
import gzip
import argparse
from liftover import ChainFile



def get_file(fname, read_mode=True):
    """get a file handler for the provided filename"""
    mode = 'r' if read_mode else 'w'
    if fname:
        if fname.endswith('.gz'):
            return gzip.open(fname, mode+'t', newline='')
        return open(fname, mode, newline='')
    return sys.stdin if read_mode else sys.stdout


def main(args):
    """read from the input, convert the CHROM and POS, then write to the output"""
    converter = ChainFile(args.chain, 'hg38', 'hg19')
    # open the files
    with get_file(args.counts) as read_file:
        with get_file(args.out, False) as write_file:
            # open the files as TSVs
            reader = csv.reader(read_file, delimiter="\t")
            writer = csv.writer(write_file, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
            # write the header verbatim, if desired
            if args.ignore_header:
                writer.writerow(next(reader))
            # iterate through each line of the input
            for line in reader:
                # convert the CHROM and POS values
                query = converter[line[0]][int(line[1])]
                # check if a liftover position exists (ie that query is not empty)
                if query:
                    assert len(query) == 1, query
                    writer.writerow(list(query[0][:2])+line[2:])


def get_args():
    """parse the arguments to this script using argparse"""
    parser = argparse.ArgumentParser(description=
        "Liftover a TSV whose first two columns are CHROM and POS, respectively."
    )
    parser.add_argument(
        '-o', '--out', default=None,
        help='output file name (in tsv format) (default: stdout)'
    )
    parser.add_argument(
        '-i', '--ignore-header', action='store_true',
        help="pass the first line (ie the header) through verbatim without liftover"
    )
    parser.add_argument(
        'chain', help="A chain file for performing the lift over"
    )
    parser.add_argument(
        'counts', default=None, nargs='?', help=
        "A TSV whose first two columns are the contig and position (default: stdin)"
    )
    args = parser.parse_args()
    return args

if __name__=="__main__":
    main(get_args())
