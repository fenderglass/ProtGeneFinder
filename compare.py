#!/usr/bin/env python

from __future__ import print_function
from find_genes import get_data_proteome, get_data_genome
import sys


E_VALUE = 0.01


#TODO: resolve problems with gene names lookup
def compare(proteome_results, genome_results, proteome_table):
    pro_rec, pro_int, pro_fam = get_data_proteome(proteome_results,
                                                  proteome_table, E_VALUE)
    gen_rec, gen_int, gen_fam = get_data_genome(genome_results, E_VALUE)

    gen_family_by_spec = {}
    for fam in gen_fam:
        for spec in fam.spectrum_ids:
            gen_family_by_spec[spec] = fam

    for f_id, p_fam in enumerate(pro_fam):
        #TODO: sort by evalue
        matched_gen_fam = None
        for pro_spec in p_fam.spectrum_ids:
            if pro_spec in gen_family_by_spec:
                matched_gen_fam = gen_family_by_spec[pro_spec]
                break

        if matched_gen_fam is not None:
            print("Matched family", f_id)
        else:
            print("Unmatched family", f_id)


def main():
    if len(sys.argv) != 4:
        print("Usage: compare.py proteome_results"
              " genome_results proteome_table", file=sys.stderr)
        return 1

    compare(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
