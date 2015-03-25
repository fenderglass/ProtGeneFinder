#!/usr/bin/env python2.7

from __future__ import print_function
import sys
from common import read_gene_matches
from collections import defaultdict
import argparse
import urllib2

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def compare_by_positions(ref_records, qry_records, only_missmatched, blast):
    #various indexes
    ref_by_fam = defaultdict(list)
    for r in ref_records:
        if r.orf_id is not None and (-1 not in [r.start, r.end]):
            ref_by_fam[r.orf_id].append(r)

    ref_family_pos = {}
    for fam, records in ref_by_fam.items():
        starts = map(lambda r: r.start, records)
        ends = map(lambda r: r.end, records)
        ref_family_pos[fam] = (min(starts), max(ends))

    qry_by_fam = defaultdict(list)
    for r in qry_records:
        if r.orf_id is not None and (-1 not in [r.start, r.end]):
            qry_by_fam[r.orf_id].append(r)

    qry_family_pos = {}
    for fam, records in qry_by_fam.items():
        starts = map(lambda r: r.start, records)
        ends = map(lambda r: r.end, records)
        qry_family_pos[fam] = (min(starts), max(ends))

    ref_rec_by_spec = {r.spec_id : r for r in ref_records}
    qry_rec_by_spec = {r.spec_id : r for r in qry_records}

    print("ORF_id\tSpec_id\tFound\tRef_pval\tRef_eval\tQry_pval\t"
          "Qry_eval\tRef_start\tQry_start\tRef_prot\tQry_prot")

    for r_fam_id, r_fam in ref_by_fam.items():
        #Matching query family based on overlap
        matches = []
        for q_fam_id in qry_by_fam:
            q_pos = qry_family_pos[q_fam_id]
            r_pos = ref_family_pos[r_fam_id]
            overlap = max(0, min(q_pos[1], r_pos[1]) - max(q_pos[0], r_pos[0]))
            if overlap > 0:
                matches.append((q_fam_id, overlap))

        matched_qry_fam = (sorted(matches, key=lambda t: t[1])[-1][0]
                           if matches else None)

        #now choose best families' spectrum
        ref_spec = min(ref_by_fam[r_fam_id], key=lambda r: r.e_value).spec_id
        if matched_qry_fam:
            qry_specs = map(lambda r: r.spec_id, qry_by_fam[matched_qry_fam])
            if ref_spec in qry_specs:
                qry_spec = ref_spec
            else:
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
            n_hits = check_blast(inner_prot)
            if n_hits:
                print("Hits: {0}".format(n_hits))
            else:
                continue

        found = "+" if matched_qry_fam is not None else "-"
        print("{0}\t{1}\t{2}\t{3:6.2e}\t{4:6.2e}\t{5}\t{6}\t"
              "{7}\t\t{8}\n{9}\t{10}\n"
              .format(r_fam_id, ref_spec, found, ref_rec.p_value,
                      ref_rec.e_value, qry_pval, qry_eval,
                      ref_start, qry_start, ref_protein, qry_protein))


def check_blast(query):
    MAX_EVAL = 0.01

    retries = 10
    records = None
    while retries > 0:
        try:
            result = NCBIWWW.qblast("blastp", "nr", query, expect=MAX_EVAL)
            records = NCBIXML.parse(result)
            break
        except urllib2.HTTPError as e:
            print(e)
            retries -= 1
    if records is None:
        raise Exception("Everything is bad!")

    n_hits = 0
    for rec in records:
        for aln in rec.alignments:
            for hsp in aln.hsps:
                n_hits += 1
    return n_hits


def main():
    parser = argparse.ArgumentParser(description="Compare two GM files",
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
    compare_by_positions(ref_gene_match, qry_gene_match, args.missmatch, args.blast)


if __name__ == "__main__":
    main()
