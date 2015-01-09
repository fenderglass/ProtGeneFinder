#!/usr/bin/env python

from __future__ import print_function
import sys
import re
from collections import namedtuple, defaultdict
from itertools import combinations
import argparse

from common import (Prsm, GeneMatch, Interval, parse_msalign_output,
                    gene_match_serialize)


def assign_intervals(records, protein_table):
    prot_table_data = {}
    with open(protein_table, "r") as f:
        for line in f:
            tokens = line.strip().split()
            start, end = int(tokens[1]), int(tokens[2])
            strand = 1 if tokens[3] == "+" else -1
            prot_table_data[tokens[0]] = (start, end, strand)

    fail_counter = 0
    for rec in records:
        prot_id = rec.prot_name.split(" ")[0].split("|")[1]
        if prot_id not in prot_table_data:
            rec.interval = Interval(-1, -1, 1)
            fail_counter += 1
            continue

        prot_rec = prot_table_data[prot_id]
        prot_len = prot_rec[1] - prot_rec[0] + 1
        assert prot_len > 0
        strand = prot_rec[2]

        if prot_rec[2] > 0:
            start = prot_rec[0] + (rec.first_res - 1) * 3
            end = prot_rec[0] + (rec.last_res - 1) * 3
        else:
            start = prot_rec[0] + (prot_len - (rec.last_res - 1) * 3 - 1)
            end = prot_rec[0] + (prot_len - (rec.first_res - 1) * 3 - 1)

        rec.interval = Interval(start, end, strand)


def assign_families(records):
    families = []
    by_prsm = {r.prsm_id : r for r in records}

    group_spectra = defaultdict(list)
    for rec in records:
        group_spectra[rec.prot_name].append(rec)

    for f_id, group in enumerate(group_spectra.values()):
        for rec in group:
            rec.family = f_id


def filter_evalue(records, e_value):
    return list(filter(lambda r: r.e_value < e_value, records))


def filter_spectras(records):
    groups = defaultdict(list)
    for rec in records:
        groups[rec.spec_id].append(rec)

    to_keep = set()
    for group in groups.itervalues():
        by_eval = sorted(group, key=lambda r: r.e_value)
        to_keep.add(by_eval[0])

    return [r for r in records if r in to_keep]


def get_matches(table_file, prot_table, e_value):
    prsms = parse_msalign_output(table_file)
    prsms = filter_spectras(prsms)
    trusted_prsms = filter_evalue(prsms, e_value)

    assign_intervals(prsms, prot_table)
    assign_families(trusted_prsms)

    matches = []
    for p in prsms:
        matches.append(GeneMatch(p.family, p.prsm_id, p.spec_id, p.p_value,
                                 p.e_value, p.interval.start,
                                 p.interval.end, p.interval.strand,
                                 p.peptide, p.genome_seq, p.html))
    return matches


def main():
    parser = argparse.ArgumentParser(description="Processing MSAlign proteome run",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("msalign_output", metavar="msalign_output",
                        help="path to result_table.txt")
    parser.add_argument("prot_table", metavar="prot_table",
                        help="path to protein table")
    parser.add_argument("-f", "--family", action="store_const",
                        dest="family", default=False, const=True,
                        help="group by families")
    parser.add_argument("-e", "--eval", dest="e_value",
                        help="custom e-value threshold",
                        default="0.01")

    args = parser.parse_args()

    gene_match = get_matches(args.msalign_output, args.prot_table,
                             float(args.e_value))
    gene_match_serialize(gene_match, sys.stdout, args.family)
    return 0


if __name__ == "__main__":
    main()
