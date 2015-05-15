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

def parse_msalign_output(filename):
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


def read_gene_matches(filename):
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


def gene_match_serialize(records, stream, family_mode):
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


def print_orf_clusters(orf_id, orf_clusters, genome_fasta, out_stream):
    #TODO: group by proteoforms
    #TODO: show modifications
    START_CODONS = ["ATG", "GTG", "TTG"]
    STOP_CODONS = ["TAG", "TAA", "TGA"]

    starts = map(lambda r: r.start, chain(*orf_clusters))
    ends = map(lambda r: r.end, chain(*orf_clusters))
    orf_begin, orf_end = min(starts), max(ends)
    sequences = get_fasta(genome_fasta)
    repr_rec = orf_clusters[0][0]
    if repr_rec.start == -1:
        return

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
    orf_seq = str(Seq(orf_seq).translate())

    def for_record(rec):
        left_flank = (rec.start - orf_begin) / 3
        right_flank = (orf_end - rec.end) / 3
        shift = ((orf_begin - left_stop - 1) / 3 if rec.strand > 0
                 else (right_stop - orf_end + 3) / 3)

        genome_seq = sequences[rec.chr_id].seq[rec.start-1:rec.end]
        if rec.strand < 0:
            genome_seq = genome_seq.reverse_complement()
            left_flank, right_flank = right_flank, left_flank
        out_stream.write("{0}{1}{2}{3}\t{4}\t{5:3.2e}\n"
                          .format(" " * shift, "." * left_flank,
                                  genome_seq.translate(), "." * right_flank,
                                  rec.spec_id, rec.e_value))

    orf_cluster_len = (orf_end - orf_begin + 1) / 3
    orf_len = (right_stop - left_stop - first_start) / 3
    strand = "+" if repr_rec.strand > 0 else "-"
    out_stream.write("ORF #{0}; Cluster length {1} aa; "
                     "ORF length {2} aa; Start {3}; End {4}; Strand {5}; "
                     "#Prsms {6}\n\n".format(orf_id, orf_cluster_len,
                     orf_len, left_stop + 3, right_stop - 3,
                     strand, len(starts)))
    out_stream.write("".join(codons_marked) + "\n")
    out_stream.write(orf_seq + "\n\n")

    for cluster in orf_clusters:
        map(lambda r: for_record(r), cluster)
        out_stream.write("\n")


def get_fasta(filename):
    if filename not in get_fasta.cache:
        get_fasta.cache[filename] = {r.id : r for r in
                                     SeqIO.parse(filename, "fasta")}
    return get_fasta.cache[filename]
get_fasta.cache = {}
