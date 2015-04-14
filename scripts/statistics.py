#!/usr/bin/env python2.7

from __future__ import print_function
import os
import sys
from collections import defaultdict, namedtuple

from Bio import SeqIO

spectrogene_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, spectrogene_root)
from spectrogene.datatypes import read_gene_matches, Interval


def calc_statistics(gene_matches, genome_file):
    sequences = {}
    for seq in SeqIO.parse(genome_file, "fasta"):
        sequences[seq.id] = seq.seq

    #all identified spectra
    print("Total spectra:\t{0}".format(len(gene_matches)))
    identified = []
    for gm in gene_matches:
        if gm.orf_id is not None:
            identified.append(Interval(gm.start, gm.end, gm.strand))

    print("\nBy identified spectra:\n------------------------")
    print("Total:\t{0} ({1:4.2f}%)"
          .format(len(identified),
                  100 * float(len(identified)) / len(gene_matches)))
    process_group(identified, sequences)
    ###

    #now group by unique peptides
    seen = set()
    uniqe_peptides = []
    for gm in identified:
        if (gm.start, gm.end) not in seen:
            seen.add((gm.start, gm.end))
            uniqe_peptides.append(Interval(gm.start, gm.end, gm.strand))

    print("\n\nBy unique proteoform\n--------------------------")
    print("Total:\t{0}".format(len(uniqe_peptides)))
    process_group(uniqe_peptides, sequences)
    ###

    #group by ORF
    by_orf = defaultdict(list)
    for gm in gene_matches:
        if gm.orf_id is not None:
            by_orf[gm.orf_id].append(gm)
    orf_intervals = []
    for orf_cluster in by_orf.values():
        start, end = sys.maxint, -sys.maxint
        for gm in orf_cluster:
            start = min(start, gm.start)
            end = max(end, gm.end)
        orf_intervals.append(Interval(start, end, orf_cluster[0].strand))
    ##########

    print("\n\nBy ORF cluster\n--------------------------")
    print("Total:\t{0}".format(len(orf_intervals)))
    process_group(orf_intervals, sequences)


def process_group(intervals, sequences):
    START_CODONS = ["ATG", "GTG", "TTG"]

    #TODO: add chr_id to GM
    chr_id = sequences.keys()[0]

    num_matched = len(intervals)
    num_start_codon_right = 0
    num_start_codon_left = 0
    num_stop_codon = 0
    num_signal = 0
    num_orf = 0
    num_stop_inside = 0
    lengths = []

    for interval in intervals:
        FLANK = 36
        if interval.start == -1:
            continue
        seq = sequences[chr_id][interval.start - 1 - FLANK :
                                interval.end + FLANK]
        if interval.strand < 0:
            seq = seq.reverse_complement()
        peptide = seq.translate()

        left = FLANK / 3
        right = left + (interval.end - interval.start) / 3

        lengths.append((interval.end - interval.start) / 3)
        start_ok = False
        stop_ok = False

        prec_codon = str(seq[FLANK -3 : FLANK])
        first_codon = str(seq[FLANK : FLANK + 3])
        if prec_codon in START_CODONS or first_codon in START_CODONS:
            num_start_codon_left += 1
            start_ok = True

        if peptide[right + 1] == "*":
            num_stop_codon += 1
            stop_ok = True

        if start_ok and stop_ok:
            num_orf += 1

        if (peptide[left - 1] == "A" and peptide[left] == "A"):
            num_signal += 1

        if "*" in peptide[left:right]:
            num_stop_inside += 1

    num_proper_start = num_start_codon_right + num_start_codon_left
    median_len = sorted(lengths)[len(lengths) / 2]

    print("Median protein len:\t{0}".format(median_len))
    print("Begins with start codon or precedes by a Start "
          "Codon (undergone NME):\t{0} ({1:4.2f}%)"
          .format(num_proper_start, 100 * float(num_proper_start) / num_matched))
    print("Ends with stop codon:\t{0} ({1:4.2f}%)"
          .format(num_stop_codon, 100 * float(num_stop_codon) / num_matched))
    print("ORF (begins/preceds with start codon AND ends "
          "with stop):\t{0} ({1:4.2f}%)"
          .format(num_orf, 100 * float(num_orf) / num_matched))
    print("Stop Codons inside: {0}".format(num_stop_inside))
    print("")
    print("Canonical signal peptide site AA:\t{0} ({1:4.2f}%)"
          .format(num_signal, 100 * float(num_signal) / num_matched))


def main():
    if len(sys.argv) != 3:
        print("Usage: statistics.py gm_file genome_file", file=sys.stderr)
        return 1
    calc_statistics(read_gene_matches(sys.argv[1]), sys.argv[2])


if __name__ == "__main__":
    main()
