from __future__ import print_function
import sys
import re
from collections import namedtuple, defaultdict
from itertools import combinations
import os
import shutil

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

    for rec in records:
        prot_id = rec.prot_name.split(" ")[0].split("|")[1]
        if prot_id not in prot_table_data:
            rec.interval = Interval(-1, -1, 1)
            continue

        p_start, p_end, p_strand = prot_table_data[prot_id]
        prot_len = p_end - p_start + 1

        if p_strand > 0:
            start = p_start + rec.first_res * 3
            end = p_start + rec.last_res * 3
        else:
            start = p_start + (prot_len - rec.last_res * 3)
            end = p_start + (prot_len - rec.first_res * 3)

        rec.interval = Interval(start, end, p_strand)


def assign_orf(records):
    families = []
    by_prsm = {r.prsm_id : r for r in records}

    group_spectra = defaultdict(list)
    for rec in records:
        group_spectra[rec.prot_name].append(rec)

    for f_id, group in enumerate(group_spectra.values()):
        for rec in group:
            rec.orf_id = f_id


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


def copy_html(prsms, out_dir):
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.mkdir(out_dir)
    for prsm in prsms:
        html_name = os.path.join(out_dir, "spec{0}.html".format(prsm.spec_id))
        shutil.copy2(prsm.html, html_name)


def get_matches(table_file, prot_table, e_value):
    prsms = parse_msalign_output(table_file)
    prsms = filter_spectras(prsms)
    trusted_prsms = filter_evalue(prsms, e_value)

    assign_intervals(prsms, prot_table)
    assign_orf(trusted_prsms)

    matches = []
    for p in prsms:
        matches.append(GeneMatch(p.orf_id, p.prsm_id, p.spec_id, p.p_value,
                                 p.e_value, p.interval.start,
                                 p.interval.end, p.interval.strand,
                                 p.peptide, p.genome_seq))
    return prsms, matches


def process_proteome(alignment_table, protein_table, e_value, out_dir):
    prsms, gene_match = get_matches(alignment_table, protein_table, e_value)

    html_dir = os.path.join(out_dir, "prsm_html")
    copy_html(prsms, html_dir)

    out_file = os.path.join(out_dir, "proteome.gm")
    gene_match_serialize(gene_match, open(out_file, "w"), False)
