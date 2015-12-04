#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

"""
This module provide data types used in SpectroGene as well as some IO functions
"""

from itertools import chain
from collections import namedtuple, defaultdict

from Bio.Seq import Seq
from Bio import SeqIO

class Prsm:
    def __init__(self, spec_id, prot_name, first_res,
                 last_res, peptide, p_value, e_value, html):
        self.spec_id = spec_id
        self.prot_name = prot_name
        self.first_res = first_res
        self.last_res = last_res
        self.peptide = peptide
        self.p_value = p_value
        self.e_value = e_value
        self.html = html

        self.orf_id = None
        self.interval = None
        self.genome_seq = None
        self.chr_id = None


Interval = namedtuple("Interval", ["start", "end", "strand"])

GeneMatch = namedtuple("GeneMatch", ["orf_id", "spec_id", "p_value", "e_value",
                                     "chr_id", "start", "end", "strand",
                                     "peptide", "genome_seq"])

def parse_toppic_output(filename):
    """
    Reads TopPic output
    """
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Data"):
                continue

            vals = line.split("\t")
            rows.append(Prsm(int(vals[2]), vals[11], int(vals[13]),
                             int(vals[14]), vals[15], float(vals[19]),
                             float(vals[20]), vals[23]))

    return rows


def parse_spectrogene_prsms(filename):
    """
    Reads SpectroGene PrSM output
    """
    gene_matches = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("ORF_id"):
                continue

            vals = line.strip().split("\t")
            strand = 1 if vals[7] == "+" else -1
            orf_id = int(vals[0]) if vals[0] != "*" else None
            genome_seq = vals[9] if vals[9] != "*" else None

            gene_matches.append(GeneMatch(orf_id, int(vals[1]), float(vals[2]),
                                float(vals[3]), vals[4], int(vals[5]),
                                int(vals[6]), strand, vals[8], genome_seq))

    return gene_matches


def write_spectrogene_prsms(records, stream, family_mode):
    """
    Writes SpectroGene PrSMs
    """
    rec_by_orf = defaultdict(list)
    without_fam = []
    for r in records:
        if r.orf_id is not None:
            rec_by_orf[r.orf_id].append(r)
        else:
            without_fam.append(r)

    stream.write("ORF_id\tSpec_id\tP_value\tE_val\tChr_id\tStart\tEnd\tStrand\t"
                 "Peptide\tGenome_seq\n")

    for orf_id in chain(rec_by_orf.values(), [without_fam]):
        by_eval = sorted(orf_id, key=lambda r: r.e_value)
        if family_mode:
            if by_eval[0].orf_id is not None:
                by_eval = [by_eval[0]]
            else:
                continue

        for m in by_eval:
            strand = "+" if m.strand > 0 else "-"
            orf_id = str(m.orf_id) if m.orf_id is not None else "*"
            genome_seq = m.genome_seq if m.genome_seq is not None else "*"

            stream.write("{0}\t{1}\t{2:4.2e}\t{3:4.2e}\t{4}\t{5}\t{6}\t{7}\t{8}"
                         "\t{9}\n"
                         .format(orf_id, m.spec_id, m.p_value, m.e_value,
                                 m.chr_id, m.start, m.end, strand, m.peptide,
                                 genome_seq))


def read_fasta(filename):
    if filename not in read_fasta.cache:
        read_fasta.cache[filename] = {r.id : r for r in
                                     SeqIO.parse(filename, "fasta")}
    return read_fasta.cache[filename]
read_fasta.cache = {}
