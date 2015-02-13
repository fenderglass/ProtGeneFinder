from itertools import chain
from collections import namedtuple, defaultdict

class Prsm:
    def __init__(self, prsm_id, spec_id, prot_name, first_res,
                 last_res, peptide, p_value, e_value, html):
        self.prsm_id = prsm_id
        self.spec_id = spec_id
        self.prot_name = prot_name
        self.first_res = first_res
        self.last_res = last_res
        self.peptide = peptide
        self.p_value = p_value
        self.e_value = e_value
        self.html = html

        self.family = None
        self.interval = None
        self.genome_seq = None


Interval = namedtuple("Interval", ["start", "end", "strand"])

GeneMatch = namedtuple("GeneMatch", ["family", "prsm_id", "spec_id", "p_value",
                                     "e_value", "start", "end", "strand",
                                     "peptide", "genome_seq", "html"])

def parse_msalign_output(filename):
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Data"):
                continue

            vals = line.split("\t")
            rows.append(Prsm(int(vals[1]), int(vals[2]), vals[11], int(vals[13]),
                             int(vals[14]), vals[15], float(vals[19]),
                             float(vals[20]), vals[23]))

    return rows


def read_gene_matches(filename):
    gene_matches = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("Fam_id"):
                continue

            vals = line.strip().split("\t")
            strand = 1 if vals[7] == "+" else -1
            family = int(vals[0]) if vals[0] != "*" else None
            genome_seq = vals[9] if vals[9] != "*" else None
            html = vals[10] if vals[10] != "*" else None

            gene_matches.append(GeneMatch(family, int(vals[1]), int(vals[2]),
                                float(vals[3]), float(vals[4]), int(vals[5]),
                                int(vals[6]), strand, vals[8], genome_seq,
                                html))

    return gene_matches


def gene_match_serialize(records, stream, family_mode):
    rec_by_fam = defaultdict(list)
    without_fam = []
    for r in records:
        if r.family is not None:
            rec_by_fam[r.family].append(r)
        else:
            without_fam.append(r)

    stream.write("Fam_id\tPrsm_id\tSpec_id\tP_value\tE_val\tStart\tEnd\tStrand\t"
                 "Peptide\tGenome_seq\tHtml\n")

    for family in chain(rec_by_fam.values(), [without_fam]):
        by_eval = sorted(family, key=lambda r: r.e_value)
        if family_mode:
            if by_eval[0].family is not None:
                by_eval = [by_eval[0]]
            else:
                continue

        for m in by_eval:
            strand = "+" if m.strand > 0 else "-"
            family = str(m.family) if m.family is not None else "*"
            genome_seq = m.genome_seq if m.genome_seq is not None else "*"

            stream.write("{0}\t{1}\t{2}\t{3:4.2e}\t{4:4.2e}\t{5}\t{6}\t{7}\t{8}"
                         "\t{9}\t{10}\n"
                         .format(family, m.prsm_id, m.spec_id, m.p_value,
                                 m.e_value, m.start, m.end, strand, m.peptide,
                                 genome_seq, m.html))


##################
class SetObject:
    pass


def MakeSet(x):
    s = SetObject()
    s.parent = s
    s.rank   = 0
    s.data = x
    return s


def Union(x, y):
    xRoot = Find(x)
    yRoot = Find(y)
    if xRoot.rank > yRoot.rank:
        yRoot.parent = xRoot
    elif xRoot.rank < yRoot.rank:
        xRoot.parent = yRoot
    elif xRoot != yRoot:
        yRoot.parent = xRoot
        xRoot.rank = xRoot.rank + 1


def Find(x):
    if x.parent == x:
       return x
    else:
       x.parent = Find(x.parent)
       return x.parent
##################
