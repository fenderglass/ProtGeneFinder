from __future__ import print_function
import sys
import re
from collections import namedtuple, defaultdict
from itertools import combinations
import os
import shutil

from Bio.Seq import Seq
from Bio import SeqIO

from common import (Prsm, GeneMatch, Interval, parse_msalign_output,
                    gene_match_serialize)


def assign_intervals(records, protein_table):
    prot_table_data = {}
    with open(protein_table, "r") as f:
        for line in f:
            tokens = line.strip().split()
            start, strand, chr_id = int(tokens[1]), tokens[2], tokens[3]
            strand = 1 if strand == "+" else -1
            prot_table_data[tokens[0]] = (start, strand, chr_id)

    for rec in records:
        prot_id = rec.prot_name.split(" ")[0].split("|")[1]
        if prot_id not in prot_table_data:
            rec.interval = Interval(-1, -1, 1)
            continue

        p_start, p_strand, p_chr_id = prot_table_data[prot_id]

        first, last = rec.first_res, rec.last_res + 1
        if p_strand > 0:
            start = p_start + first * 3
            end = p_start + last * 3 - 1
        else:
            start = p_start - last * 3 + 1
            end = p_start - first * 3

        rec.interval = Interval(start, end, p_strand)
        rec.chr_id = p_chr_id
        rec.prot_id = prot_id


"""
def assign_genome_seqs(records, genome_file):
    FLANK_LEN = 10
    genome = get_fasta(genome_file)

    for record in records:
        if record.interval.start == -1:
            continue

        flank_start = (record.interval.start - 1) - (FLANK_LEN * 3)
        flank_end = record.interval.end + (FLANK_LEN * 3)
        genome_seq = genome[record.chr_id].seq[flank_start:flank_end]

        if record.interval.strand < 0:
            genome_seq = genome_seq.reverse_complement()

        translated = str(genome_seq.translate())
        translated = ".".join([translated[0:FLANK_LEN].lower(),
                               translated[FLANK_LEN:-FLANK_LEN],
                               translated[-FLANK_LEN:].lower()])

        record.genome_seq = translated


def get_fasta(filename):
    return {r.id : r for r in SeqIO.parse(filename, "fasta")}
"""


def assign_orf(records):
    families = []

    group_spectra = defaultdict(list)
    for rec in records:
        group_spectra[rec.prot_name].append(rec)

    for f_id, group in enumerate(group_spectra.values()):
        for rec in group:
            rec.orf_id = f_id + 1


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
    #assign_genome_seqs(prsms, "datasets/Salmonella_msalign/S.Typhimurium.fasta")
    assign_orf(trusted_prsms)

    matches = []
    for p in prsms:
        matches.append(GeneMatch(p.orf_id, p.spec_id, p.p_value, p.e_value,
                                 p.chr_id, p.interval.start, p.interval.end,
                                 p.interval.strand, p.peptide, p.genome_seq))
    return prsms, matches


def process_proteome(alignment_table, protein_table, e_value, out_dir):
    prsms, gene_match = get_matches(alignment_table, protein_table, e_value)

    html_dir = os.path.join(out_dir, "prsm_html")
    copy_html(prsms, html_dir)

    out_file = os.path.join(out_dir, "proteome.gm")
    gene_match_serialize(gene_match, open(out_file, "w"), False)
