#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

"""
This module outputs ORF clusters into a text file
"""

from collections import defaultdict

from Bio.Seq import Seq

from spectrogene.datatypes import read_fasta


def print_orf_clusters(orf_id, prsms, genome_file, out_stream):
    """
    Pretty ORF cluster printer
    """
    START_CODONS = ["ATG", "GTG", "TTG"]
    STOP_CODONS = ["TAG", "TAA", "TGA"]
    CLEAV_AA = "GASCTPV"
    NO_CLEAV_AA = "DEFHIKLMNOQRUWY"

    starts = map(lambda r: r.start, prsms)
    ends = map(lambda r: r.end, prsms)
    orf_begin, orf_end = min(starts), max(ends)
    sequences = read_fasta(genome_file)
    repr_rec = prsms[0]
    if repr_rec.start == -1:
        return

    #group by proteoforms
    by_proteoform = defaultdict(list)
    for rec in prsms:
        by_proteoform[(rec.start, rec.end)].append(rec)
    proteoforms = []
    for pf_list in by_proteoform.values():
        proteoforms.append(min(pf_list, key=lambda r: r.e_value))
    ##

    #getting flanking stop codons
    def get_orf_flank(seq, orf_begin, strand):
        left_stop = orf_begin - 3
        while True:
            codon = seq[left_stop: left_stop + 3]
            if strand < 0:
                codon = str(Seq(codon).reverse_complement())
            if codon in STOP_CODONS or left_stop < 0:
                break
            left_stop -= 3
        right_stop = orf_begin + 3
        while True:
            codon = chr_seq[right_stop: right_stop + 3]
            if strand < 0:
                codon = str(Seq(codon).reverse_complement())
            if codon in STOP_CODONS or right_stop >= len(chr_seq) - 3:
                break
            right_stop += 3
        return left_stop, right_stop

    ###getting ORF sequence and marking potential start codons
    chr_seq = str(sequences[repr_rec.chr_id].seq)
    left_stop, right_stop = get_orf_flank(chr_seq, orf_begin - 1,
                                          repr_rec.strand)
    orf_seq = chr_seq[left_stop : right_stop + 3]
    if repr_rec.strand < 0:
        orf_seq = str(Seq(orf_seq).reverse_complement())
    codons_marked = ""
    first_start = None
    for i in xrange(len(orf_seq) / 3):
        codon = orf_seq[i * 3:(i + 1) * 3]
        codons_marked += "v" if codon in START_CODONS else " "
        if codon in START_CODONS and not first_start:
            first_start = i * 3
    orf_translated = str(Seq(orf_seq).translate())
    #######

    #some statistics
    cluster_shift = (orf_begin - left_stop - 1 if repr_rec.strand > 0
                     else right_stop - orf_end + 3)
    beg_start = orf_seq[cluster_shift : cluster_shift + 3] in START_CODONS
    prec_start = orf_seq[cluster_shift - 3 : cluster_shift] in START_CODONS
    nme_good = False if beg_start or prec_start else None
    if beg_start and orf_translated[cluster_shift / 3 + 1] in NO_CLEAV_AA:
        nme_good = True
    if prec_start and orf_translated[cluster_shift / 3] in CLEAV_AA:
        nme_good = True
    weak_signal = (cluster_shift / 3 >= 1 and
                   orf_translated[cluster_shift / 3 - 1] == "A")
    strong_signal = (weak_signal and cluster_shift / 3 >= 3 and
                     orf_translated[cluster_shift / 3 - 3] == "A")

    orf_cluster_len = (orf_end - orf_begin + 1) / 3
    orf_len = (right_stop - left_stop - first_start) / 3
    strand = repr_rec.strand > 0
    ###

    out_stream.write("=" * len(codons_marked) + "\n\n")
    out_stream.write("ORF#: {0}\nCluster length: {1} aa\n"
                     "ORF length: {2} aa\n#Prsms: {6}\n"
                     "#Proteoforms: {10}\n"
                     "Start: {3} nt\nEnd: {4} nt\nStrand: {5}\n"
                     "Begins with Start: {7}\nPreceeds with Start: {8}\n"
                     "NME: {9}\nSignal peptide motif (weak): {11}\n"
                     "Signal peptide morif (strong): {12}\n\n"
                     .format(orf_id, orf_cluster_len, orf_len,
                     left_stop + 3, right_stop - 3, _bs(strand), len(prsms),
                     _bs(beg_start), _bs(prec_start), _bs(nme_good),
                      len(proteoforms), _bs(weak_signal), _bs(strong_signal)))
    out_stream.write("".join(codons_marked) + "\n")
    out_stream.write(orf_translated + "\n\n")

    for rec in sorted(proteoforms, key=lambda p: p.e_value):
        #for each record
        left_flank = (rec.start - orf_begin) / 3
        right_flank = (orf_end - rec.end) / 3
        genome_seq = sequences[rec.chr_id].seq[rec.start-1:rec.end]
        mod_str = _modification_string(rec)
        if rec.strand < 0:
            genome_seq = genome_seq.reverse_complement()
            left_flank, right_flank = right_flank, left_flank
        out_stream.write("{0}{1}\n".format(" " * ((cluster_shift / 3) + left_flank),
                                           mod_str))
        out_stream.write("{0}{1}{2}{3}\t{4}\t{5:3.2e}\n"
                          .format(" " * (cluster_shift / 3), "." * left_flank,
                                  genome_seq.translate(), "." * right_flank,
                                  rec.spec_id, rec.e_value))
    out_stream.write("\n")


def _modification_string(prsm_rec):
    """
    Constructs a string for PTM annotation
    """
    prsm_len = (prsm_rec.end - prsm_rec.start) / 3
    dot_l, dot_r = prsm_rec.peptide.find("."), prsm_rec.peptide.rfind(".")
    mod_str = ""
    cur_mod = False
    mod_weight = False
    cur_weight = ""
    shift = 0
    for ch in prsm_rec.peptide[dot_l + 1 : dot_r]:
        if shift != 0:
            shift -= 1
        elif ch == "(":
            cur_mod = True
        elif ch == ")":
            cur_mod = False
        elif ch == "[":
            mod_weight = True
        elif ch == "]":
            weight = ("42" if cur_weight == "Acetylation" else
                      str(int(float(cur_weight))))
            mod_str += weight
            shift = len(weight)
            mod_weight = False
            cur_weight = ""
        elif cur_mod:
            mod_str += "_"
        elif mod_weight:
            cur_weight += ch
        else:
            mod_str += " "

    return mod_str


def _bs(x):
    """
    Returns the sign of the argument
    """
    if x is None:
        return "NA"
    return "+" if x else "-"
