#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

"""
This module creates ORFeome
"""

from __future__ import print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


MIN_ORF_LEN = 10


def _chunks(string, size):
    """
    Splits the string into chunks of the given size
    """
    for i in range(0, len(string), size):
        yield i, string[i:i+size]


def _orf_partition(sequence):
    """
    Partitions sequence into ORFs (by Stop codons)
    """
    START_CODONS = ["ATG", "GTG", "TTG"]
    STOP_CODONS = ["TAG", "TAA", "TGA"]

    cur_orf = 0
    sequence += STOP_CODONS[0]
    for i in xrange(0, len(sequence), 3):
        codon = str(sequence[i : i + 3])
        if codon not in STOP_CODONS:
            continue

        for start_cand in xrange(cur_orf, i, 3):
            codon = str(sequence[start_cand : start_cand + 3])
            if codon in START_CODONS:
                yield start_cand, sequence[start_cand : i].translate()
                break

        cur_orf = i + 3


def make_proteome(filename, window_size, out_proteome):
    """
    The main function
    """
    out_fasta = open(out_proteome, "w")

    for record in SeqIO.parse(filename, "fasta"):
        #forward
        for frame_shift in range(0, 3):
            trim = (len(record.seq) - frame_shift) % 3
            genome_seq = (record.seq[frame_shift:-trim] if trim else
                          record.seq[frame_shift:])

            #partition into ORFs
            for orf_pos, orf in _orf_partition(genome_seq):
                for ovlp in [0, window_size / 2]:
                    for win_pos, window in _chunks(orf[ovlp:], window_size):
                        if len(window) < MIN_ORF_LEN:
                            continue

                        shift = 3 * (win_pos + ovlp) + orf_pos + frame_shift
                        window_name = "{0}::fwd_{1}".format(record.id, shift)
                        SeqIO.write(SeqRecord(seq=window, id=window_name,
                                    description=""), out_fasta, "fasta")

        #reverse
        for frame_shift in range(0, 3):
            trim = (len(record.seq) - frame_shift) % 3
            rev_seq = record.seq.reverse_complement()
            genome_seq = (rev_seq[frame_shift:-trim] if trim else
                          rev_seq[frame_shift:])

            #partition into ORFs
            for orf_pos, orf in _orf_partition(genome_seq):
                for ovlp in [0, window_size / 2]:
                    for win_pos, window in _chunks(orf[ovlp:], window_size):
                        if len(window) < MIN_ORF_LEN:
                            continue

                        shift = 3 * (win_pos + ovlp) + orf_pos + frame_shift
                        shift = len(record.seq) - 1 - shift
                        window_name = "{0}::rev_{1}".format(record.id, shift)
                        SeqIO.write(SeqRecord(seq=window, id=window_name,
                                    description=""), out_fasta, "fasta")
