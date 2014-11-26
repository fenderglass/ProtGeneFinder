#!/usr/bin/env python

from __future__ import print_function
import sys
from common import read_gene_matches
from collections import defaultdict
import argparse


#TODO: resolve problems with gene names lookup
def compare(ref_records, qry_records, only_missmatched):
    ref_by_fam = defaultdict(list)
    for r in ref_records:
        if r.family is not None:
            ref_by_fam[r.family].append(r)

    qry_by_fam = defaultdict(list)
    for r in ref_records:
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
        pro_matches = sorted(r_fam, key=lambda m: m.e_value)

        for pro_match in pro_matches:
            if qry_rec_by_spec[pro_match.spec_id].family is not None:
                matched_qry_fam = qry_rec_by_spec[pro_match.spec_id].family
                break

        pro_spec = pro_match.spec_id
        if matched_qry_fam and only_missmatched:
            continue

        ref_rec = ref_rec_by_spec[pro_spec]
        qry_rec = qry_rec_by_spec[pro_spec]

        start_ref = str(ref_rec_by_spec[pro_spec].start)
        start_qry = str(qry_rec_by_spec[pro_spec].start)
        found = "+" if matched_qry_fam else "-"

        start_ref = start_ref if start_ref != "-1" else "n/a"
        start_qry = start_qry if start_qry != "-1" else "n/a"

        ref_protein = (ref_rec.genome_seq if ref_rec.genome_seq else
                       ref_rec.peptide)
        qry_protein = (qry_rec.genome_seq if qry_rec.genome_seq else
                       qry_rec.peptide)

        print("{0}\t{1}\t{2}\t{3:6.2e}\t{4:6.2e}\t{5:6.2e}\t{6:6.2e}\t"
              "{7}\t\t{8}\n{9}\t{10}"
              .format(r_fam_id, pro_spec, found, ref_rec.p_value,
                      ref_rec.e_value, qry_rec.p_value, qry_rec.e_value,
                      start_ref, start_qry, ref_protein, qry_protein))

        print(ref_rec.html)
        print("")


def main():
    parser = argparse.ArgumentParser(description="Compare two sets of genome matches",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("ref_table", metavar="ref_table",
                        help="path to reference table in gm fromat")
    parser.add_argument("qry_table", metavar="qry_table",
                        help="path to querry table in gm format")
    parser.add_argument("-m", "--missmatch", action="store_const",
                        dest="missmatch", default=False, const=True,
                        help="show only missmatches")

    args = parser.parse_args()

    ref_gene_match = read_gene_matches(args.ref_table)
    qry_gene_match = read_gene_matches(args.qry_table)
    compare(ref_gene_match, qry_gene_match, args.missmatch)


if __name__ == "__main__":
    main()
