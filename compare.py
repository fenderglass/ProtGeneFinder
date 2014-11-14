#!/usr/bin/env python

from __future__ import print_function
from find_genes import get_data_proteome, get_data_genome
import sys


E_VALUE = 0.01
PROT_REF = True


#TODO: resolve problems with gene names lookup
def compare(ref_records, ref_families, qry_records, qry_families):

    ###SOME MAGIC INDEXES
    #TODO: check if one spectrum points to different families (like repeat)
    ref_rec_by_prsm = {r.prsm_id : r for r in ref_records}
    ref_rec_by_spec = {}
    for f in ref_families:
        for prsm in f.prsms:
            rec = ref_rec_by_prsm[prsm]
            ref_rec_by_spec[rec.spec_id] = rec

    #TODO: check if one spectrum points to different families (like repeat)
    qry_rec_by_prsm = {r.prsm_id : r for r in qry_records}
    qry_family_by_spec = {}
    qry_rec_by_spec = {}
    for fam in qry_families:
        for prsm in fam.prsms:
            rec = qry_rec_by_prsm[prsm]
            qry_family_by_spec[rec.spec_id] = fam
            qry_rec_by_spec[rec.spec_id] = rec
    #add spectrums that are not in families (low e-value)
    for rec in qry_records:
        if rec.spec_id not in qry_rec_by_spec:
            qry_rec_by_spec[rec.spec_id] = rec
    ######################

    print("Fam_id\tFound\tRef_pval\tRef_eval\tQry_pval\t"
          "Qry_eval\tRef_start\tQry_start")
    for f_num, r_fam in enumerate(ref_families):
        matched_qry_fam = None
        pro_prsms = sorted(r_fam.prsms, key=lambda p: ref_rec_by_prsm[p].e_value)
        pro_specs = [ref_rec_by_prsm[p].spec_id for p in pro_prsms]

        for pro_spec in pro_specs:
            if pro_spec in qry_family_by_spec:
                matched_qry_fam = qry_family_by_spec[pro_spec]
                break

        ref_rec = ref_rec_by_spec[pro_spec]
        qry_rec = qry_rec_by_spec[pro_spec]
        #p_val_ref = ref_rec_by_spec[pro_spec].p_value
        #e_val_ref = ref_rec_by_spec[pro_spec].e_value
        #p_val_qry = qry_rec_by_spec[pro_spec].p_value
        #e_val_qry = qry_rec_by_spec[pro_spec].e_value

        start_ref = str(r_fam.start)
        if matched_qry_fam is None:
            start_qry = str(qry_rec_by_spec[pro_spec].interval.start)
            found = "-"
        else:
            start_qry = str(matched_qry_fam.start)
            found = "+"

        start_ref = start_ref if start_ref != "-1" else "n/a"
        start_qry = start_qry if start_qry != "-1" else "n/a"

        print("{0}\t{1}\t{2:6.2e}\t{3:6.2e}\t{4:6.2e}\t{5:6.2e}\t{6}\t\t{7}"
              .format(f_num, found, ref_rec.p_value, ref_rec.e_value,
                      qry_rec.p_value, qry_rec.e_value, start_ref, start_qry))


def main():
    if len(sys.argv) != 4:
        print("Usage: compare.py proteome_results"
              " genome_results proteome_table", file=sys.stderr)
        return 1

    proteome_results = sys.argv[1]
    genome_results = sys.argv[2]
    proteome_table = sys.argv[3]

    pro_rec, pro_fam = get_data_proteome(proteome_results, proteome_table,
                                         E_VALUE)
    gen_rec, gen_fam = get_data_genome(genome_results, E_VALUE)

    if PROT_REF:
        compare(pro_rec, pro_fam, gen_rec, gen_fam)
    else:
        compare(gen_rec, gen_fam, pro_rec, pro_fam)


if __name__ == "__main__":
    main()
