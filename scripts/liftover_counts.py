#!/usr/bin/env python
import sys
import csv
from liftover import ChainFile

converter = ChainFile(sys.argv[2], 'hg38', 'hg19')

with open(sys.argv[1], 'w', newline='') as read_file:
    with open(sys.stdout, 'w', newline='') as write_file:
        writer = csv.writer(write_file, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        for line in csv.reader(csvfile, delimiter="\t"):
            print(converter[line[0], line[1]])
            # writer.write_row(line[:1])
