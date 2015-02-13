#!/usr/bin/env python

from __future__ import print_function
import sys
from common import read_gene_matches
from collections import defaultdict
import argparse
import subprocess


def compare_by_positions(ref_records, qry_records, only_missmatched, blast):
    #various indexes
    ref_by_fam = defaultdict(list)
    for r in ref_records:
        if r.family is not None and (-1 not in [r.start, r.end]):
            ref_by_fam[r.family].append(r)

    ref_family_pos = {}
    for fam, records in ref_by_fam.items():
        starts = map(lambda r: r.start, records)
        ends = map(lambda r: r.end, records)
        ref_family_pos[fam] = (min(starts), max(ends))

    qry_by_fam = defaultdict(list)
    for r in qry_records:
        if r.family is not None and (-1 not in [r.start, r.end]):
            qry_by_fam[r.family].append(r)

    qry_family_pos = {}
    for fam, records in qry_by_fam.items():
        starts = map(lambda r: r.start, records)
        ends = map(lambda r: r.end, records)
        qry_family_pos[fam] = (min(starts), max(ends))

    ref_rec_by_spec = {r.spec_id : r for r in ref_records}
    qry_rec_by_spec = {r.spec_id : r for r in qry_records}

    print("Fam_id\tSpec_id\tFound\tRef_pval\tRef_eval\tQry_pval\t"
          "Qry_eval\tRef_start\tQry_start\tRef_prot\tQry_prot")

    #matched_families = set()
    for r_fam_id, r_fam in ref_by_fam.items():
        #Matching query family based on overlap
        matches = []
        for q_fam_id in qry_by_fam:
            #if q_fam_id in matched_families: continue
            q_pos = qry_family_pos[q_fam_id]
            r_pos = ref_family_pos[r_fam_id]
            overlap = max(0, min(q_pos[1], r_pos[1]) - max(q_pos[0], r_pos[0]))
            if overlap > 0:
                matches.append((q_fam_id, overlap))

        matched_qry_fam = (sorted(matches, key=lambda t: t[1])[-1][0]
                           if matches else None)
        #matched_families.add(matched_qry_fam)

        #now choose best families' spectrum
        ref_spec = min(ref_by_fam[r_fam_id], key=lambda r: r.e_value).spec_id
        if matched_qry_fam:
            qry_spec = min(qry_by_fam[matched_qry_fam],
                           key=lambda r: r.e_value).spec_id
        else:
            qry_spec = ref_spec

        if matched_qry_fam is not None and only_missmatched:
            continue

        ref_rec = ref_rec_by_spec[ref_spec]
        qry_rec = qry_rec_by_spec.get(qry_spec, None)

        ref_start = str(ref_rec.start)
        ref_protein = (ref_rec.genome_seq if ref_rec.genome_seq else
                       ref_rec.peptide)

        if qry_rec:
            qry_start = str(qry_rec.start)
            qry_protein = (qry_rec.genome_seq if qry_rec.genome_seq else
                           qry_rec.peptide)
            qry_pval = "%6.2e" % qry_rec.p_value
            qry_eval = "%6.2e" % qry_rec.e_value
        else:
            qry_start, qry_protein, qry_pval, qry_eval = "n/a", "n/a", "n/a", "n/a"

        if blast:
            left = ref_protein.find(".")
            right = ref_protein.rfind(".")
            inner_prot = ref_protein[left+1:right]
            if only_missmatched and not check_blast(inner_prot):
                continue

        found = "+" if matched_qry_fam is not None else "-"
        print("{0}\t{1}\t{2}\t{3:6.2e}\t{4:6.2e}\t{5}\t{6}\t"
              "{7}\t\t{8}\n{9}\t{10}"
              .format(r_fam_id, ref_spec, found, ref_rec.p_value,
                      ref_rec.e_value, qry_pval, qry_eval,
                      ref_start, qry_start, ref_protein, qry_protein))

        print("")


def compare_by_spectrum(ref_records, qry_records, only_missmatched, blast):
    ref_by_fam = defaultdict(list)
    for r in ref_records:
        if r.family is not None:
            ref_by_fam[r.family].append(r)

    qry_by_fam = defaultdict(list)
    for r in qry_records:
        if r.family is not None:
            qry_by_fam[r.family].append(r)

    ref_rec_by_prsm = {r.prsm_id : r for r in ref_records}
    ref_rec_by_spec = {r.spec_id : r for r in ref_records}

    qry_rec_by_prsm = {r.prsm_id : r for r in qry_records}
    qry_rec_by_spec = {r.spec_id : r for r in qry_records}

    print("Fam_id\tSpec_id\tFound\tRef_pval\tRef_eval\tQry_pval\t"
          "Qry_eval\tRef_start\tQry_start\tRef_prot\tQry_prot")
    for r_fam_id, r_fam in ref_by_fam.items():
        matched_qry_fam = None
        ref_matches = sorted(r_fam, key=lambda m: m.e_value)

        for ref_match in ref_matches:
            if qry_rec_by_spec[ref_match.spec_id].family is not None:
                matched_qry_fam = qry_rec_by_spec[ref_match.spec_id].family
                break

        ref_spec = ref_match.spec_id
        qry_spec = ref_spec
        if matched_qry_fam is not None and only_missmatched:
            continue

        ref_rec = ref_rec_by_spec[ref_spec]
        qry_rec = qry_rec_by_spec[qry_spec]

        start_ref = str(ref_rec_by_spec[ref_spec].start)
        start_qry = str(qry_rec_by_spec[qry_spec].start)
        found = "+" if matched_qry_fam is not None else "-"

        start_ref = start_ref if start_ref != "-1" else "n/a"
        start_qry = start_qry if start_qry != "-1" else "n/a"

        ref_protein = (ref_rec.genome_seq if ref_rec.genome_seq else
                       ref_rec.peptide)
        qry_protein = (qry_rec.genome_seq if qry_rec.genome_seq else
                       qry_rec.peptide)

        if blast:
            left = ref_protein.find(".")
            right = ref_protein.rfind(".")
            inner_prot = ref_protein[left+1:right]
            if only_missmatched and not check_blast(inner_prot):
                continue

        print("{0}\t{1}\t{2}\t{3:6.2e}\t{4:6.2e}\t{5:6.2e}\t{6:6.2e}\t"
              "{7}\t\t{8}\n{9}\t{10}"
              .format(r_fam_id, ref_spec, found, ref_rec.p_value,
                      ref_rec.e_value, qry_rec.p_value, qry_rec.e_value,
                      start_ref, start_qry, ref_protein, qry_protein))

        print("")


def check_blast(query):
    MAX_EVAL = 0.01

    proc = subprocess.Popen(["blastp", "-db", "nr", "-remote", "-outfmt",
                             "6 sseqid evalue", "-max_target_seqs", "1"],
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE)
    stdout, _stderr = proc.communicate(query)
    if not stdout or float(stdout.split("\n")[0].split("\t")[1]) > MAX_EVAL:
        return False
    return True


def _median(values):
    """
    Not a true median, but we keep real distances
    """
    sorted_values = sorted(values)
    return sorted_values[(len(values) - 1) / 2]


def main():
    parser = argparse.ArgumentParser(description="Compare two sets of genome matches",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("ref_table", metavar="ref_table",
                        help="path to reference table in gm format")
    parser.add_argument("qry_table", metavar="qry_table",
                        help="path to querry table in gm format")
    parser.add_argument("-m", "--missmatch", action="store_const",
                        dest="missmatch", default=False, const=True,
                        help="show only missmatches")
    parser.add_argument("-b", "--blast", action="store_const",
                        dest="blast", default=False, const=True,
                        help="check missmatches with blast")

    args = parser.parse_args()

    ref_gene_match = read_gene_matches(args.ref_table)
    qry_gene_match = read_gene_matches(args.qry_table)
    #compare_by_spectrum(ref_gene_match, qry_gene_match, args.missmatch, args.blast)
    compare_by_positions(ref_gene_match, qry_gene_match, args.missmatch, args.blast)


if __name__ == "__main__":
    main()
