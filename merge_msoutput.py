#!/usr/bin/env python

import sys

def merge_tables(tables):
    header = False
    prsm_id = 0
    for table in tables:
        with open(table, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("Data"):
                    if not header:
                        print(line)
                        header = True
                    continue

                vals = line.split("\t")
                vals[1] = str(prsm_id)
                prsm_id += 1
                print("\t".join(vals))


def main():
    merge_tables(sys.argv[1:])


if __name__ == "__main__":
    main()
