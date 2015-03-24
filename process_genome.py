#!/usr/bin/env python

from __future__ import print_function
import sys
import re
from collections import defaultdict
from itertools import combinations
import shutil
import os

from Bio.Seq import Seq
from Bio import SeqIO

from common import (Prsm, GeneMatch, Interval, parse_msalign_output,
                    gene_match_serialize, MakeSet, Find, Union)


def assign_intervals(records):
    CONV_SHIFT = 1
    for rec in records:
        seq_name, meta = rec.prot_name.split(" ")[0].split("::")
        direction, shift_len, genome_pos = meta.split("_")

        rec.seq_name = seq_name

        if direction == "fwd":
            genomic_start = int(genome_pos) + rec.first_res * 3
            genomic_end = int(genome_pos) + rec.last_res * 3
        elif direction == "rev":
            genomic_start = int(genome_pos) - rec.last_res * 3 + 1    #why +1??
            genomic_end = int(genome_pos) - rec.first_res * 3 + 1     #why +1??

        assert genomic_end >= genomic_start

        rec.interval = Interval(genomic_start + CONV_SHIFT,
                                genomic_end + CONV_SHIFT,
                                1 if direction == "fwd" else -1)


def assign_genome_seqs(records, genome_file):
    FLANK_LEN = 20
    genome = get_fasta(genome_file)

    for record in records:
        if record.interval is None:
            continue

        flank_start = (record.interval.start - 1) - (FLANK_LEN * 3)
        flank_end = (record.interval.end - 1) + (FLANK_LEN * 3)
        genome_seq = genome[record.seq_name].seq[flank_start:flank_end]

        if record.interval.strand < 0:
            genome_seq = genome_seq.reverse_complement()

        genome_seq += "N" * (-len(genome_seq) % 3)

        translated = str(genome_seq.translate())
        translated = ".".join([translated[0:FLANK_LEN].lower(),
                               translated[FLANK_LEN:-FLANK_LEN+1],
                               translated[-FLANK_LEN+1:].lower()])
        record.genome_seq = translated


def assign_families(records, genome_file):
    #genome = get_fasta(genome_file)

    sets = {r.prsm_id : MakeSet(r) for r in records}
    for rec_1, rec_2 in combinations(records, 2):
        int_1 = rec_1.interval
        int_2 = rec_2.interval
        if rec_1.seq_name != rec_2.seq_name:
            continue

        #test for overlapping
        overlap = (min(int_1.end, int_2.end) -
                   max(int_1.start, int_2.start))
        if overlap > 0:
            Union(sets[rec_1.prsm_id], sets[rec_2.prsm_id])
        """
        #"linking"
        elif abs(overlap) < 90 and abs(int_1.start - int_2.start) % 3 == 0:
            gap_start = min(int_1.end, int_2.end)
            gap_end = max(int_1.start, int_2.start)
            gap_seq = genome[rec_1.seq_name].seq[gap_start:gap_end].translate()
            if "*" not in gap_seq:
                Union(sets[rec_1.prsm_id], sets[rec_2.prsm_id])
        """

    by_family = defaultdict(list)
    for s in sets.values():
        by_family[Find(s)].append(s.data)

    for fam_id, prsms in enumerate(by_family.values()):
        for prsm in prsms:
            prsm.family = fam_id


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


def get_fasta(filename):
    return {r.id : r for r in SeqIO.parse(filename, "fasta")}


def copy_html(prsms, out_dir):
    for prsm in prsms:
        html_name = os.path.join(out_dir, "spec{0}.html".format(prsm.spec_id))
        shutil.copy2(prsm.html, html_name)


def get_matches(table_file, genome_file, e_value):
    prsms = parse_msalign_output(table_file)
    prsms = filter_spectras(prsms)
    trusted_prsms = filter_evalue(prsms, e_value)


    assign_intervals(prsms)
    assign_genome_seqs(prsms, genome_file)
    assign_families(trusted_prsms, genome_file)

    matches = []
    for p in prsms:
        matches.append(GeneMatch(p.family, p.prsm_id, p.spec_id, p.p_value,
                                 p.e_value, p.interval.start,
                                 p.interval.end, p.interval.strand,
                                 p.peptide, p.genome_seq))

    return prsms, matches


def process_genome(alignment_table, fasta_file, evalue, out_dir):
    prsms, gene_match = get_matches(alignment_table, fasta_file,
                                    evalue)

    html_dir = os.path.join(out_dir, "prsm_html")
    if os.path.isdir(html_dir):
        shutil.rmtree(html_dir)
    os.mkdir(html_dir)
    copy_html(prsms, html_dir)

    out_file = os.path.join(out_dir, "genome.gm")
    gene_match_serialize(gene_match, open(out_file, "w"), False)
