from common_processor import CommonProcessor
from common import Interval

class ProteomeProcessor(CommonProcessor):
    def __init__(self, e_value, genome_fasta, prot_alignment):
        super(ProteomeProcessor, self).__init__(e_value, genome_fasta)
        self.prot_alignment = prot_alignment

    def assign_intervals(self):
        prot_table_data = {}
        with open(self.prot_alignment, "r") as f:
            for line in f:
                tokens = line.strip().split()
                start, strand, chr_id = int(tokens[1]), tokens[2], tokens[3]
                strand = 1 if strand == "+" else -1
                prot_table_data[tokens[0]] = (start, strand, chr_id)

        for rec in self.prsms:
            prot_id = rec.prot_name.split(" ")[0].split("|")[1]
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
