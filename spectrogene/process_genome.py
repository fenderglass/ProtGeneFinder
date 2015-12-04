#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

"""
This class extends CommonProcessor for the genome run case
"""

from spectrogene.common_processor import CommonProcessor
from spectrogene.datatypes import Interval

class GenomeProcessor(CommonProcessor):
    def __init__(self, e_value, genome_fasta):
        super(GenomeProcessor, self).__init__(e_value, genome_fasta)

    def _assign_intervals(self):
        """
        A function for assigning genomic coordinate to PrSMs
        """
        CONV_SHIFT = 1
        for rec in self.prsms:
            seq_name, meta = rec.prot_name.split(" ")[0].split("::")
            direction, genome_pos = meta.split("_")

            rec.chr_id = seq_name

            first, last = rec.first_res, rec.last_res + 1
            if direction == "fwd":
                genomic_start = int(genome_pos) + first * 3 + CONV_SHIFT
                genomic_end = int(genome_pos) + last * 3

            elif direction == "rev":
                genomic_start = int(genome_pos) - last * 3 + 2
                genomic_end = int(genome_pos) - first * 3 + 1

            rec.interval = Interval(genomic_start, genomic_end,
                                    1 if direction == "fwd" else -1)
