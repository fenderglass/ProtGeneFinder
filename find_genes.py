#!/usr/bin/env python

import sys
from collections import namedtuple, defaultdict
from itertools import combinations


E_VALUE = 0.01

Prsm = namedtuple("Prsm", ["spec_id", "prot_name", "first_res", "last_res",
                           "peptide", "e_value"])
Interval = namedtuple("Interval", ["spectrum_id", "start", "end", "strand"])
Family = namedtuple("Family", ["spectrum_ids", "start", "end"])

def parse_table(filename):
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Data"):
                continue

            vals = line.split()
            rows.append(Prsm(int(vals[2]), vals[9], int(vals[11]),
                             int(vals[12]), vals[13], float(vals[18])))

    return rows


def get_intervals(records):
    intervals = []
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

        intervals.append(Interval(rec.spec_id, genomic_start, genomic_end,
                                  1 if direction == "fwd" else -1))

    return intervals


def filter_evalue(records, e_value):
    return list(filter(lambda r: r.e_value < e_value, records))


def get_families(intervals):
    sets = {i.spectrum_id : MakeSet(i.spectrum_id) for i in intervals}
    for int_1, int_2 in combinations(intervals, 2):
        #test for overlapping
        if ((int_1.start <= int_2.start and int_2.end <= int_1.end) or
            (int_2.start <= int_1.start and int_1.end <= int_2.end)):
            Union(sets[int_1.spectrum_id], sets[int_2.spectrum_id])

    by_family = defaultdict(list)
    for s in sets.values():
        by_family[Find(s)].append(s.data)
    by_spectrum = {}
    for i in intervals:
        by_spectrum[i.spectrum_id] = i

    families = []
    for fam in by_family.values():
        spectrums = fam
        start = min([by_spectrum[i].start for i in spectrums])
        end = max([by_spectrum[i].end for i in spectrums])
        families.append(Family(spectrums, start, end))

    return families


def print_table(records, intervals, families):
    rec_by_id = {rec.spec_id : rec for rec in records}
    int_by_id = {i.spectrum_id : i for i in intervals}
    for i, family in enumerate(sorted(families, key=lambda f: f.start)):
        prot_len = int((family.end - family.start) / 3)
        print("Family {0}\t{1}\t{2}".format(i, family.start, prot_len))

        by_eval = sorted(family.spectrum_ids, key=lambda p: rec_by_id[p].e_value)
        for spec in by_eval:
            interval = int_by_id[spec]
            record = rec_by_id[spec]
            prot_len = int((interval.end - interval.start) / 3)
            sign = "+" if interval.strand > 0 else "-"
            print("\t{0}\t{1}\t{2}\t{3:4.2e}\t{4}\t{5}"
                    .format(interval.start, prot_len, sign,
                            record.e_value, record.spec_id, record.peptide))

        print("")


def get_data(table_files, e_value):
    records = []
    for file in table_files:
        records.extend(parse_table(file))

    records = filter_evalue(records, e_value)
    intervals = get_intervals(records)
    families = get_families(intervals)
    return records, intervals, families


def main():
    table_files = sys.argv[1:]
    print_table(*get_data(table_files, E_VALUE))


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
