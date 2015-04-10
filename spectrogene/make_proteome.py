from __future__ import print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def chunks(string, size):
    for i in range(0, len(string), size):
        yield i, string[i:i+size]


def make_proteome(filename, window_size, out_proteome):
    MIN_ORF = 10
    out_fasta = open(out_proteome, "w")

    for record in SeqIO.parse(filename, "fasta"):
        #forward
        for frame_shift in range(0, 3):
            trim = (len(record.seq) - frame_shift) % 3
            prot_seq = (record.seq[frame_shift:-trim].translate() if trim else
                        record.seq[frame_shift:].translate())

            #partition into ORFs
            orf_shift = frame_shift
            frames = str(prot_seq).split("*")
            for orf in frames:
                for ovlp in [0, window_size / 2]:
                    for win_pos, window in chunks(orf[ovlp:], window_size):
                        if len(window) < MIN_ORF:
                            continue

                        shift = 3 * (win_pos + ovlp) + orf_shift
                        window_name = "{0}::fwd_{1}".format(record.id, shift)
                        SeqIO.write(SeqRecord(seq=Seq(window), id=window_name,
                                    description=""), out_fasta, "fasta")

                orf_shift += (len(orf) + 1) * 3

        #reverse
        for frame_shift in range(0, 3):
            trim = (len(record.seq) - frame_shift) % 3
            rev_seq = record.seq.reverse_complement()
            prot_seq = (rev_seq[frame_shift:-trim].translate() if trim else
                        rev_seq[frame_shift:].translate())

            #partition into ORFs
            orf_shift = frame_shift
            frames = str(prot_seq).split("*")
            for orf in frames:
                for ovlp in [0, window_size / 2]:
                    for win_pos, window in chunks(orf[ovlp:], window_size):
                        if len(window) < MIN_ORF:
                            continue

                        shift = 3 * (win_pos + ovlp) + orf_shift
                        shift = len(record.seq) - 1 - shift
                        window_name = "{0}::rev_{1}".format(record.id, shift)
                        SeqIO.write(SeqRecord(seq=Seq(window), id=window_name,
                                    description=""), out_fasta, "fasta")

                orf_shift += (len(orf) + 1) * 3
