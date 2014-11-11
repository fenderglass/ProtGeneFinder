#!/usr/bin/env python

from __future__ import print_function
import sys

def make_table(ncbi_table, uniprot_table):
    prot_table_data = {}
    with open(ncbi_table, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            vals = line.strip().split("\t")
            start, end = int(vals[2]), int(vals[3])
            strand = 1 if vals[4] == "+" else -1
            stm = vals[7]

            prot_table_data[stm] = (start, end, strand)

    with open(uniprot_table, "r") as f:
        for line in f:
            if line.startswith("Entry"):
                continue
            vals = line.strip().split("\t")
            if " " in vals[4]:
                gene_id = vals[4].split(" ")[0]
                stm = vals[4].split(" ")[-1]
                if len(vals[4].split(" ")) > 2:
                    print (vals[4])
                #print(gene_id, *prot_table_data[stm])
            #else:
                #print(vals[4], *prot_table_data[vals[4]])


def main():
    if len(sys.argv) != 3:
        print("Usage: make_prototable.py ncbi_table uniprot_table")
        return 1

    make_table(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
