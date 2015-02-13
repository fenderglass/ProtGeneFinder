#!/usr/bin/env python

from __future__ import print_function
import os
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


PROT_FILENAME_1 = "proteome_1.fasta"
PROT_FILENAME_2 = "proteome_2.fasta"


def chunks(string, size):
    for i in range(0, len(string), size):
        yield i, string[i:i+size]


def make_proteome(filename, slice_size, out_dir):
    assert slice_size % 2 == 0
    non_shift = os.path.join(out_dir, PROT_FILENAME_1)
    slice_prot(filename, slice_size, 0, open(non_shift, "w"))

    half_shift = os.path.join(out_dir, PROT_FILENAME_2)
    slice_prot(filename, slice_size, int(slice_size / 2), open(half_shift, "w"))

    return os.path.abspath(non_shift), os.path.abspath(half_shift)


def slice_prot(filename, slice_size, glob_shift, file_out):
    for record in SeqIO.parse(filename, "fasta"):
        #forward
        for frame_shift in range(0, 3):
            shift = frame_shift + glob_shift * 3
            trim = (len(record.seq) - shift) % 3
            prot_seq = (record.seq[shift:-trim].translate() if trim else
                        record.seq[shift:].translate())

            for prot_pos, chunk in chunks(prot_seq, slice_size):
                genome_pos = 3 * prot_pos + shift
                prot_name = "{0}::fwd_{1}_{2}".format(record.id, str(frame_shift),
                                                      genome_pos)
                chunk = Seq(str(chunk).replace("*", "G"))
                SeqIO.write(SeqRecord(seq=chunk, id=prot_name, description=""),
                            file_out, "fasta")

        #reverse
        for frame_shift in range(0, 3):
            shift = frame_shift + glob_shift * 3
            trim = (len(record.seq) - shift) % 3
            prot_seq = (record.seq.reverse_complement()[shift:-trim].translate()
                        if trim else
                        record.seq.reverse_complement()[shift:].translate())

            for prot_pos, chunk in chunks(prot_seq, slice_size):
                genome_pos = len(record.seq) - 1 - (3 * prot_pos + shift)
                prot_name = "{0}::rev_{1}_{2}".format(record.id, str(frame_shift),
                                                      genome_pos)
                chunk = Seq(str(chunk).replace("*", "X"))
                SeqIO.write(SeqRecord(seq=chunk, id=prot_name, description=""),
                            file_out, "fasta")


def main():
    if len(sys.argv) != 4:
        print("Usage: make_proteome.py dna_fasta slice_size out_dir",
              file=sys.stderr)
        return 1

    slice_size = int(sys.argv[2])
    if slice_size % 2 == 1:
        print("Slice size must be even", file=sys.stderr)
        return 1

    make_proteome(sys.argv[1], int(sys.argv[2]), sys.argv[3])
    return 0


if __name__ == "__main__":
    main()
