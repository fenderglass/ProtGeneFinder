#!/usr/bin/env python

from __future__ import print_function
from find_genes import get_data_proteome, get_data_genome
import sys


E_VALUE = 0.01


#TODO: resolve problems with gene names lookup
def compare(proteome_results, genome_results, proteome_table):
    pro_rec, pro_fam = get_data_proteome(proteome_results, proteome_table,
                                         E_VALUE)
    gen_rec, gen_fam = get_data_genome(genome_results, E_VALUE)

    ###SOME MAGIC INDEXES
    #TODO: check if one spectrum points to different families (like repeat)
    pro_rec_by_prsm = {r.prsm_id : r for r in pro_rec}
    pro_rec_by_spec = {}
    for f in pro_fam:
        for prsm in f.prsms:
            rec = pro_rec_by_prsm[prsm]
            pro_rec_by_spec[rec.spec_id] = rec

    #TODO: check if one spectrum points to different families (like repeat)
    gen_rec_by_prsm = {r.prsm_id : r for r in gen_rec}
    gen_family_by_spec = {}
    gen_rec_by_spec = {}
    for fam in gen_fam:
        for prsm in fam.prsms:
            rec = gen_rec_by_prsm[prsm]
            gen_family_by_spec[rec.spec_id] = fam
            gen_rec_by_spec[rec.spec_id] = rec
    #add spectrums that are not in families (low e-value)
    for rec in gen_rec:
        if rec.spec_id not in gen_rec_by_spec:
            gen_rec_by_spec[rec.spec_id] = rec
    ######################

    print("Fam_id\tFound\tProt_pval\tGen_pval\tProt_start\tGen_start")
    for f_num, p_fam in enumerate(pro_fam):
        matched_gen_fam = None
        pro_prsms = sorted(p_fam.prsms, key=lambda p: pro_rec_by_prsm[p].e_value)
        pro_specs = [pro_rec_by_prsm[p].spec_id for p in pro_prsms]

        for pro_spec in pro_specs:
            if pro_spec in gen_family_by_spec:
                matched_gen_fam = gen_family_by_spec[pro_spec]
                break

        if matched_gen_fam is None:
            p_val_pro = pro_rec_by_spec[pro_spec].p_value
            p_val_gen = gen_rec_by_spec[pro_spec].p_value
            start_pro = str(p_fam.start) if p_fam.start != -1 else "n/a"
            start_gen = gen_rec_by_spec[pro_spec].interval.start
            print("{0}\t-\t{1:6.2e}\t{2:6.2e}\t{3}\t\t{4}"
                  .format(f_num, p_val_pro, p_val_gen, start_pro, start_gen))
        else:
            p_val_pro = pro_rec_by_spec[pro_spec].p_value
            p_val_gen = gen_rec_by_spec[pro_spec].p_value
            start_pro = str(p_fam.start) if p_fam.start != -1 else "n/a"
            start_gen = matched_gen_fam.start
            print("{0}\t+\t{1:6.2e}\t{2:6.2e}\t{3}\t\t{4}"
                  .format(f_num, p_val_pro, p_val_gen, start_pro, start_gen))


def main():
    if len(sys.argv) != 4:
        print("Usage: compare.py proteome_results"
              " genome_results proteome_table", file=sys.stderr)
        return 1

    compare(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
