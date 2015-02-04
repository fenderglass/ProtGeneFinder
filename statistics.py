#!/usr/bin/env python

from __future__ import print_function
import sys
from collections import defaultdict

from common import read_gene_matches

def calc_statistics(gene_matches):
    print("Total spectra:\t{0}".format(len(gene_matches)))
    ident_spectra = [g for g in gene_matches if g.family is not None]

    print("\nBy identified spectra:\n------------------------")
    print("Total:\t{0} ({1:4.2f}%)"
          .format(len(ident_spectra),
                  100 * float(len(ident_spectra)) / len(gene_matches)))
    process_group(ident_spectra)

    #now group by unique peptides
    seen = set()
    uniqe_peptides = []
    for gm in ident_spectra:
        if (gm.start, gm.end) not in seen:
            seen.add((gm.start, gm.end))
            uniqe_peptides.append(gm)

    print("\n\nBy unique identified peptide\n--------------------------")
    print("Total:\t{0}".format(len(uniqe_peptides)))
    process_group(uniqe_peptides)

    #group by family
    by_family = defaultdict(list)
    for gm in ident_spectra:
        by_family[gm.family].append(gm)
    longest_in_fam = []
    for fam in by_family.values():
        longest_in_fam.append(sorted(fam, key=lambda g: g.e_value)[0])

    print("\n\nBy family (for best match)\n--------------------------")
    print("Total:\t{0}".format(len(longest_in_fam)))
    process_group(longest_in_fam)


def process_group(gene_matches):
    #START_CODONS = ["M", "V", "L"]
    START_CODONS = ["M"]
    num_matched = len(gene_matches)
    num_start_codon_right = 0
    num_start_codon_left = 0
    num_stop_codon = 0
    num_signal = 0
    num_orf = 0
    num_suspicious = 0
    lengths = []
    stop_inside = 0
    stop_orf = 0

    for match in gene_matches:
        peptide = match.genome_seq if match.genome_seq else match.peptide
        left = peptide.find(".")
        right = peptide.rfind(".")

        lengths.append(right - left)
        start_ok = False
        stop_ok = False

        if left > 0 and peptide[left - 1].upper() in START_CODONS:
            num_start_codon_left += 1
            start_ok = True
        if peptide[left + 1].upper() in START_CODONS:
            num_start_codon_right += 1
            start_ok = True

        if len(peptide) > right + 1 and peptide[right + 1] == "*":
            num_stop_codon += 1
            stop_ok = True

        if start_ok and stop_ok:
            if "*" in peptide[left + 1:right]:
                stop_orf += 1
            num_orf += 1

        if not start_ok or not stop_ok:
            start_sites = [i for i, l in enumerate(peptide) if l in START_CODONS]
            stop_sites = [i for i, l in enumerate(peptide) if l == "*"]
            if start_sites and stop_sites and min(start_sites) < max(stop_sites):
                num_suspicious += 1
        #if not (start_ok or stop_ok):
        #    if "*" in peptide[right + 1:].upper():
        #        for st in START_CODONS:
        #            if st in peptide[:left+1].upper():
        #                num_suspicious += 1
        #                break

        if (left > 0 and peptide[left - 1].upper() == "A" and
            peptide[left + 1].upper() == "A"):
            num_signal += 1

        if "*" in peptide[left + 1:right]:
            stop_inside += 1

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
    #print("extended ORF (start/stop codons "
    #      "beyond PrSM in <= 20 aa):\t{0} ({1:4.2f}%)"
    print("ORF inside PrSM:\t{0} ({1:4.2f}%)"
          .format(num_suspicious, 100 * float(num_suspicious) / num_matched))
    print("")
    print("Stop codon inside protein:\t{0} ({1:4.2f}%)"
          .format(stop_inside, 100 * float(stop_inside) / num_matched))
    print("Stop codon inside AND has exact ORF:\t{0} ({1:4.2f}%)"
          .format(stop_orf, 100 * float(stop_orf) / num_matched))
    print("")
    print("Canonical signal peptide site AA:\t{0} ({1:4.2f}%)"
          .format(num_signal, 100 * float(num_signal) / num_matched))


def main():
    if len(sys.argv) != 2:
        print("Usage: statistics.py gm_file", file=sys.stderr)
    calc_statistics(read_gene_matches(sys.argv[1]))


if __name__ == "__main__":
    main()
