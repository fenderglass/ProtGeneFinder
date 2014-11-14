#!/usr/bin/env python

from __future__ import print_function
import sys
import re
from collections import namedtuple, defaultdict
from itertools import combinations


E_VALUE = 0.01

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

#Prsm = namedtuple("Prsm", ["spec_id", "prot_name", "first_res", "last_res",
#                           "peptide", "p_value", "e_value"])
Interval = namedtuple("Interval", ["start", "end", "strand"])
Family = namedtuple("Family", ["prsms", "start", "end"])

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


def get_intervals_genome(records):
    CONV_SHIFT = 1
    for rec in records:
        meta = rec.prot_name.split("::")[1]
        direction, shift_len, genome_pos = meta.split("_")

        if direction == "fwd":
            genomic_start = int(genome_pos) + (rec.first_res - 1) * 3
            genomic_end = int(genome_pos) + (rec.last_res - 1) * 3
        elif direction == "rev":
            genomic_start = int(genome_pos) - (rec.last_res - 1) * 3
            genomic_end = int(genome_pos) - (rec.first_res - 1) * 3

        assert genomic_end >= genomic_start

        rec.interval = Interval(genomic_start + CONV_SHIFT,
                                genomic_end + CONV_SHIFT,
                                1 if direction == "fwd" else -1)


def get_intervals_proteome(records, protein_table):
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

        start = prot_table_data[prot_id][0]
        end = prot_table_data[prot_id][1]
        strand = prot_table_data[prot_id][2]
        rec.interval = Interval(start, end, strand)

    print("Total fails:", fail_counter)


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
        families.append(Family(prsm_ids, start, end))

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
        families.append(Family(prsms, start, end))

    return families


def filter_evalue(records, e_value):
    return list(filter(lambda r: r.e_value < e_value, records))


def print_table(records, families):
    rec_by_prsm = {rec.prsm_id : rec for rec in records}
    for i, family in enumerate(sorted(families, key=lambda f: f.start)):
        prot_len = int((family.end - family.start) / 3)
        print("Family {0}\t{1}\t{2}".format(i, family.start, prot_len))

        by_eval = sorted(family.prsms, key=lambda p: rec_by_prsm[p].e_value)
        for prsm in by_eval:
            record = rec_by_prsm[prsm]
            prot_len = int((record.interval.end - record.interval.start) / 3)
            sign = "+" if record.interval.strand > 0 else "-"
            print("\t{0}\t{1}\t{2}\t{3:4.2e}\t{4}\t{5}"
                    .format(record.interval.start, prot_len, sign,
                            record.e_value, record.spec_id, record.peptide))

        print("")


def get_data_genome(table_file, e_value):
    records = parse_table(table_file)
    get_intervals_genome(records)

    filtered_records = filter_evalue(records, e_value)
    families = get_families_genome(filtered_records)

    return records, families


def get_data_proteome(results_file, prot_table, e_value):
    records = parse_table(results_file)
    get_intervals_proteome(records, prot_table)

    filtered_records = filter_evalue(records, e_value)
    families = get_families_proteome(filtered_records)
    return records, families


#TODO: collapse similar spectras
def main():
    print_table(*get_data_genome(sys.argv[1], E_VALUE))


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
