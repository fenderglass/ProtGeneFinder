#!/usr/bin/env python

from find_genes import get_data
import sys


E_VALUE = 0.01


def compare(master_table, slave_table):
    ms_records, ms_intervals, ms_families = get_data(master_table, E_VALUE)
    sl_records, sl_intervals, sl_families = get_data(slave_table, E_VALUE)

    sl_family_by_spec = {}
    for fam in sl_families:
        for spec in fam.spectrum_ids:
            sl_family_by_spec[spec] = fam

    for fam_id, ms_fam in enumerate(ms_families):
        #TODO: sort by evalue
        matched_slave_fam = None
        for ms_spec in ms_fam.spectrum_ids:
            if ms_spec in sl_family_by_spec:
                matched_slave_fam = sl_family_by_spec[ms_spec]
                break

        if matched_slave_fam is not None:
            print("Matched family", fam_id)
        else:
            print("Unmatched family", fam_id)


def main():
    compare(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
