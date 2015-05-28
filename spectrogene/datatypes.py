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


def print_orf_clusters(orf_id, prsms, genome_fasta, out_stream):
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
    sequences = get_fasta(genome_fasta)
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
    weak_signal = orf_translated[cluster_shift / 3 - 1] == "A"
    strong_signal = weak_signal and orf_translated[cluster_shift / 3 - 3] == "A"

    orf_cluster_len = (orf_end - orf_begin + 1) / 3
    orf_len = (right_stop - left_stop - first_start) / 3
    strand = repr_rec.strand > 0
    ###

    ###True gene starts
    #true_start = None
    #table = _parse_blast_alignment("datasets/Salmonella_proteome/aln.txt")
    #for (start, end, _strand) in table.values():
    #    if left_stop <= start and end <= right_stop:
    #        if _strand > 0:
    #            true_start = (start - left_stop) / 3
    #        else:
    #            true_start = (right_stop - end) / 3 + 1
    #        break
    #if true_start is not None:
    #    codons_marked = codons_marked[:true_start] + "$" + codons_marked[true_start + 1:]
    ###

    #if true_start in [cluster_shift / 3, cluster_shift / 3 - 1]:
    #    return
    #if strong_signal:
    #    return
    #if not (beg_start or prec_start):
    #    return

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

"""
def _parse_blast_alignment(filename):
    #temp
    table = {}
    with open(filename, "r") as f:
        processsed = set()
        for line in f:
            tokens = line.strip().split()
            qry_name, qry_chr, qry_start, qry_end, ref_start, ref_end  = \
                    (tokens[0], tokens[1], int(tokens[6]), int(tokens[7]),
                     int(tokens[8]), int(tokens[9]))
            if qry_name in processsed:
                continue

            processsed.add(qry_name)
            strand = 1 if ref_start < ref_end else -1
            shift = -(qry_start - 1) * 3
            if strand < 0:
                shift = -shift

            seq_id = qry_name.split("|")[1]
            length = abs(qry_start - qry_end) * 3
            if strand > 0:
                table[seq_id] = (ref_start + shift, ref_start + shift + length, strand)
            else:
                table[seq_id] = (ref_start + shift - length, ref_start + shift, strand)

    return table
"""


def _bs(x):
    if x is None:
        return "NA"
    return "+" if x else "-"


def get_fasta(filename):
    if filename not in get_fasta.cache:
        get_fasta.cache[filename] = {r.id : r for r in
                                     SeqIO.parse(filename, "fasta")}
    return get_fasta.cache[filename]
get_fasta.cache = {}
