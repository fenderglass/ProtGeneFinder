#!/usr/bin/env python

import sys
import os

def merge_tables(tables, out_stream):
    header = False
    prsm_id = 0
    for table in tables:
        directory = os.path.abspath(os.path.dirname(table))
        with open(table, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("Data"):
                    if not header:
                        out_stream.write(line + "\tHtml\n")
                        header = True
                    continue

                vals = line.split("\t")
                old_prsm_id = vals[1]
                vals[1] = str(prsm_id)
                prsm_id += 1

                html = os.path.join(directory, "spectra_html", "prsms",
                                    "prsm{0}.html".format(old_prsm_id))
                vals.append(html)

                out_stream.write("\t".join(vals) + "\n")


def main():
    merge_tables(sys.argv[1:], sys.stdout)


if __name__ == "__main__":
    main()
