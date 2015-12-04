#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

from spectrogene.common_processor import CommonProcessor
from spectrogene.datatypes import Interval

class ProteomeProcessor(CommonProcessor):
    def __init__(self, e_value, genome_fasta, prot_alignment):
        super(ProteomeProcessor, self).__init__(e_value, genome_fasta)
        self.prot_alignment = prot_alignment

    def assign_intervals(self):
        prot_table_data = _parse_blast_alignment(self.prot_alignment)

        for rec in self.prsms:
            prot_id = rec.prot_name.split(" ")[0]
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


def _parse_blast_alignment(filename):
    table = {}
    with open(filename, "r") as f:
        processsed = set()
        for line in f:
            tokens = line.strip().split()
            qry_name, qry_chr, qry_start, ref_start, ref_end  = \
                    (tokens[0], tokens[1], int(tokens[6]),
                     int(tokens[8]), int(tokens[9]))
            if qry_name in processsed:
                continue

            processsed.add(qry_name)
            strand = 1 if ref_start < ref_end else -1
            shift = -(qry_start - 1) * 3
            if strand < 0:
                shift = -shift

            #seq_id = qry_name.split("|")[1]
            table[qry_name] = (ref_start + shift, strand, qry_chr)

    return table
