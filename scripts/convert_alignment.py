#!/usr/bin/env python

from __future__ import print_function
import sys

with open(sys.argv[1], "r") as f:
    processsed = set()
    for line in f:
        tokens = line.strip().split()
        qry_name, qry_chr, qry_start, qry_end  = (tokens[0], tokens[1],
                                                 int(tokens[8]), int(tokens[9]))
        if qry_name in processsed:
            continue

        processsed.add(qry_name)
        strand = "+" if qry_start < qry_end else "-"
        if strand == "-":
            qry_start, qry_end = qry_end, qry_start
        print(qry_name.split("|")[1], qry_start, qry_end, strand, qry_chr)
