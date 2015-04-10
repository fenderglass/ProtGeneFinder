from common_processor import CommonProcessor
from datatypes import Interval

class GenomeProcessor(CommonProcessor):
    def __init__(self, e_value, genome_fasta):
        super(GenomeProcessor, self).__init__(e_value, genome_fasta)

    def assign_intervals(self):
        CONV_SHIFT = 1
        for rec in self.prsms:
            seq_name, meta = rec.prot_name.split(" ")[0].split("::")
            direction, genome_pos = meta.split("_")

            rec.chr_id = seq_name

            first, last = rec.first_res, rec.last_res + 1
            if direction == "fwd":
                genomic_start = int(genome_pos) + first * 3 + CONV_SHIFT
                genomic_end = int(genome_pos) + last * 3
                #assert (genomic_end - genomic_start + 1) % 3 == 0
            elif direction == "rev":
                genomic_start = int(genome_pos) - last * 3 + 2
                genomic_end = int(genome_pos) - first * 3 + 1
                #assert (genomic_end - genomic_start + 1) % 3 == 0

            rec.interval = Interval(genomic_start, genomic_end,
                                    1 if direction == "fwd" else -1)
