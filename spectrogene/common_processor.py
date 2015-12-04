#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

"""
This module represents an abstract class
with main SpectroGene functions
"""

from __future__ import print_function
from collections import namedtuple, defaultdict
from itertools import combinations
import os
import shutil

from spectrogene.datatypes import (GeneMatch, parse_msalign_output,
                                   gene_match_serialize, get_fasta,
                                   print_orf_clusters)
from spectrogene.disjoint_set import MakeSet, Union, Find


class CommonProcessor(object):
    def __init__(self, e_value, genome_fasta):
        self.e_value = e_value
        self.genome_fasta = genome_fasta
        self.prsms = None
        self.matches = None

    def process_and_output(self, toppic_table, out_file):
        prsms = parse_msalign_output(toppic_table)
        prsms = _keep_best_spectra(prsms)
        self.prsms = prsms

        self.assign_intervals()
        self.assign_orf()
        self.assign_genome_seqs()

        matches = []
        for p in prsms:
            matches.append(GeneMatch(p.orf_id, p.spec_id, p.p_value, p.e_value,
                                     p.chr_id, p.interval.start, p.interval.end,
                                     p.interval.strand, p.peptide, p.genome_seq))

        gene_match_serialize(matches, open(out_file, "w"), False)
        self.matches = matches

    def assign_intervals(self):
        pass

    def assign_genome_seqs(self):
        FLANK_LEN = 10
        sequences = get_fasta(self.genome_fasta)

        for record in self.prsms:
            if record.interval.start == -1:
                continue

            flank_start = (record.interval.start - 1) - (FLANK_LEN * 3)
            flank_end = record.interval.end + (FLANK_LEN * 3)
            genome_seq = sequences[record.chr_id].seq[flank_start:flank_end]

            if record.interval.strand < 0:
                genome_seq = genome_seq.reverse_complement()

            translated = str(genome_seq.translate())
            translated = ".".join([translated[0:FLANK_LEN].lower(),
                                   translated[FLANK_LEN:-FLANK_LEN],
                                   translated[-FLANK_LEN:].lower()])

            record.genome_seq = translated

    def assign_orf(self):
        MAX_GAP = 3000  #arbitrary
        sequences = get_fasta(self.genome_fasta)

        trusted_prsms = _filter_evalue(self.prsms, self.e_value)
        sets = {r : MakeSet(r) for r in trusted_prsms}

        for rec_1, rec_2 in combinations(trusted_prsms, 2):
            int_1 = rec_1.interval
            int_2 = rec_2.interval
            if rec_1.chr_id is None or rec_1.chr_id != rec_2.chr_id:
                continue

            #linking into ORFs
            overlap = (min(int_1.end, int_2.end) -
                       max(int_1.start, int_2.start) + 1)
            if overlap > 0 and overlap % 3 == 0:
                Union(sets[rec_1], sets[rec_2])

            #TODO: optimize search
            elif (abs(overlap) < MAX_GAP and
                  abs(int_1.start - int_2.start) % 3 == 0):
                gap_start = min(int_1.end, int_2.end)
                gap_end = max(int_1.start, int_2.start) - 1
                gap_seq = sequences[rec_1.chr_id].seq[gap_start:gap_end]

                if int_1.strand < 0:
                    gap_seq = gap_seq.reverse_complement()

                if "*" not in gap_seq.translate():
                    Union(sets[rec_1], sets[rec_2])

        by_orf = defaultdict(list)
        for s in sets.values():
            by_orf[Find(s)].append(s.data)

        for orf_id, prsms in enumerate(by_orf.values()):
            for prsm in prsms:
                prsm.orf_id = orf_id + 1

    def print_orfs(self, out_file):
        """
        Pretty printing of the ORF structure
        """
        out_stream = open(out_file, "w")
        by_orf = defaultdict(list)
        proteoform_len = {}
        for r in self.matches:
            if r.orf_id is not None:
                by_orf[r.orf_id].append(r)
                length = r.end - r.start
                proteoform_len[r.orf_id] = max(proteoform_len.get(r.orf_id, 0),
                                               length)

        for orf_id in sorted(by_orf, key=proteoform_len.get, reverse=True):
            print_orf_clusters(orf_id, by_orf[orf_id], self.genome_fasta,
                               out_stream)


def _filter_evalue(prsms, e_value):
    return list(filter(lambda r: r.e_value < e_value, prsms))


def _keep_best_spectra(prsms):
    groups = defaultdict(list)
    for rec in prsms:
        groups[rec.spec_id].append(rec)

    to_keep = set()
    for group in groups.itervalues():
        by_eval = sorted(group, key=lambda r: r.p_value)
        to_keep.add(by_eval[0])

    return [r for r in prsms if r in to_keep]
