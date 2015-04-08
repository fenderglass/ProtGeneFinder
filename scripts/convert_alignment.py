#!/usr/bin/env python

from __future__ import print_function
import sys

with open(sys.argv[1], "r") as f:
    processsed = set()
    for line in f:
        tokens = line.strip().split()
        qry_name, qry_chr, qry_start, ref_start, ref_end  = \
                (tokens[0], tokens[1], int(tokens[6]),
                 int(tokens[8]), int(tokens[9]))
        if qry_name in processsed:
            continue

        processsed.add(qry_name)
        strand = "+" if ref_start < ref_end else "-"
        shift = -(qry_start - 1) * 3
        if strand == "-":
            shift = -shift
        print(qry_name.split("|")[1], ref_start + shift,
              strand, qry_chr)
