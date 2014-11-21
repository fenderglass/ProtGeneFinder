#!/usr/bin/env python

from __future__ import print_function
import sys
import re
from collections import namedtuple, defaultdict
from itertools import combinations

from Bio.Seq import Seq
from Bio import SeqIO

E_VALUE = 0.01
ONLY_BEST = True

class Prsm:
    def __init__(self, prsm_id, spec_id, prot_name, first_res,
                 last_res, peptide, p_value, e_value):
        self.prsm_id = prsm_id
        self.spec_id = spec_id
        self.prot_name = prot_name
        self.first_res = first_res
        self.last_res = last_res
        self.peptide = peptide
        self.p_value = p_value
        self.e_value = e_value
        self.interval = None

Interval = namedtuple("Interval", ["start", "end", "strand"])
Family = namedtuple("Family", ["id", "prsms", "start", "end"])

def parse_table(filename):
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Data"):
                continue

            vals = line.split("\t")
            rows.append(Prsm(int(vals[1]), int(vals[2]), vals[9], int(vals[11]),
                             int(vals[12]), vals[13], float(vals[17]),
                             float(vals[18])))

    return rows


def assign_intervals_genome(records):
    CONV_SHIFT = 1
    for rec in records:
        meta = rec.prot_name.split("::")[1]
        direction, shift_len, genome_pos = meta.split("_")

        if direction == "fwd":
            genomic_start = int(genome_pos) + (rec.first_res - 1) * 3
            genomic_end = int(genome_pos) + (rec.last_res - 1) * 3
        elif direction == "rev":
            genomic_start = int(genome_pos) - (rec.last_res - 1) * 3 + 1    #why +1??
            genomic_end = int(genome_pos) - (rec.first_res - 1) * 3 + 1     #why +1??

        assert genomic_end >= genomic_start

        rec.interval = Interval(genomic_start + CONV_SHIFT,
                                genomic_end + CONV_SHIFT,
                                1 if direction == "fwd" else -1)


def assign_intervals_proteome(records, protein_table):
    prot_table_data = {}
    with open(protein_table, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            vals = line.strip().split("\t")
            start, end = int(vals[2]), int(vals[3])
            strand = 1 if vals[4] == "+" else -1
            if vals[6] != "-":
                prot_id = vals[6]
            else:
                prot_id = vals[7].split(".")[0] #mwahaha

            prot_table_data[prot_id] = (start, end, strand)

    fail_counter = 0
    gen_id_re = re.compile(".*(^|\s)GN=(\S*)($|\s).*")
    for rec in records:
        prot_id = gen_id_re.match(rec.prot_name).group(2)
        if prot_id not in prot_table_data:
            rec.interval = Interval(-1, -1, 1)
            fail_counter += 1
            continue

        prot_rec = prot_table_data[prot_id]
        prot_len = prot_rec[1] - prot_rec[0] + 1
        assert prot_len > 0
        strand = prot_rec[2]

        if prot_rec[2] > 0:
            start = prot_rec[0] + (rec.first_res - 1) * 3
            end = prot_rec[0] + (rec.last_res - 1) * 3
        else:
            start = prot_rec[0] + (prot_len - (rec.last_res - 1) * 3 - 1)
            end = prot_rec[0] + (prot_len - (rec.first_res - 1) * 3 - 1)

        rec.interval = Interval(start, end, strand)

    #print("Total fails:", fail_counter)


def get_families_proteome(records):
    families = []
    by_prsm = {r.prsm_id : r for r in records}

    group_spectra = defaultdict(list)
    for rec in records:
        group_spectra[rec.prot_name].append(rec)

    for group in group_spectra.values():
        start = sys.maxint
        end = 0
        prsm_ids = list(map(lambda r: r.prsm_id, group))
        for rec in group:
            start = min(start, by_prsm[rec.prsm_id].interval.start)
            end = max(end, by_prsm[rec.prsm_id].interval.end)
        families.append(Family(start, prsm_ids, start, end))

    return families


def get_families_genome(records):
    sets = {r.prsm_id : MakeSet(r.prsm_id) for r in records}
    for rec_1, rec_2 in combinations(records, 2):
        int_1 = rec_1.interval
        int_2 = rec_2.interval
        #test for overlapping
        if ((int_1.start <= int_2.start and int_2.end <= int_1.end) or
            (int_2.start <= int_1.start and int_1.end <= int_2.end)):
            Union(sets[rec_1.prsm_id], sets[rec_2.prsm_id])

    by_family = defaultdict(list)
    for s in sets.values():
        by_family[Find(s)].append(s.data)
    by_prsm = {r.prsm_id : r for r in records}

    families = []
    for prsms in by_family.values():
        start = min([by_prsm[p].interval.start for p in prsms])
        end = max([by_prsm[p].interval.end for p in prsms])
        families.append(Family(start, prsms, start, end))

    return families


def filter_evalue(records, e_value):
    return list(filter(lambda r: r.e_value < e_value, records))


def print_table(records, families, genome_file, only_best):
    FLANK_LEN = 5
    rec_by_prsm = {rec.prsm_id : rec for rec in records}
    genome = get_genome(genome_file)

    print("Fam_id\tSpec_id\tE_value\t\tStart\tEnd\tStrand\tPeptide")
    for family in sorted(families, key=lambda f: f.id):
        by_eval = sorted(family.prsms, key=lambda p: rec_by_prsm[p].e_value)
        if only_best:
            by_eval = [by_eval[0]]

        for prsm in by_eval:
            record = rec_by_prsm[prsm]
            strand = "+" if record.interval.strand > 0 else "-"

            ##
            seq_name = record.prot_name.split("::")[0]
            flank_start = (record.interval.start - 1) - (FLANK_LEN * 3)
            flank_end = (record.interval.end - 1) + (FLANK_LEN * 3)
            genome_seq = genome[seq_name].seq[flank_start:flank_end]

            if strand == "-":
                genome_seq = genome_seq.reverse_complement()
            translated = str(genome_seq.translate())
            translated = ".".join([translated[0:FLANK_LEN],
                                   translated[FLANK_LEN:-FLANK_LEN+1],
                                   translated[-FLANK_LEN+1:]])
            ##

            print("{0}\t{1}\t{2:4.2e}\t{3}\t{4}\t{5}\t{6}"
                    .format(family.id, record.spec_id, record.e_value,
                            record.interval.start, record.interval.end,
                            strand, translated))



def get_data_genome(table_file, e_value):
    records = parse_table(table_file)
    records = filter_spectras(records)
    assign_intervals_genome(records)

    filtered_records = filter_evalue(records, e_value)
    families = get_families_genome(filtered_records)

    return records, families


def get_data_proteome(results_file, prot_table, e_value):
    records = parse_table(results_file)
    assign_intervals_proteome(records, prot_table)

    filtered_records = filter_evalue(records, e_value)
    families = get_families_proteome(filtered_records)
    return records, families


def filter_spectras(records):
    groups = defaultdict(list)
    for rec in records:
        groups[rec.spec_id].append(rec)

    to_keep = set()
    for group in groups.itervalues():
        by_eval = sorted(group, key=lambda r: r.e_value)
        to_keep.add(by_eval[0])

    return [r for r in records if r in to_keep]


def get_genome(filename):
    return {r.id : r for r in SeqIO.parse(filename, "fasta")}


def main():
    if len(sys.argv) != 3:
        print("Usage: find_genes.py results_table genome_file")
        return 1

    records, families = get_data_genome(sys.argv[1], E_VALUE)
    print_table(records, families, sys.argv[2], ONLY_BEST)
    return 0


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


if __name__ == "__main__":
    main()
