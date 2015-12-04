#!/usr/bin/env python

#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

"""
This script merges multipe files with top-down spectra in
msalign fromat into one
"""

from __future__ import print_function
import sys

def merge(files_list):
    counter = 0
    for file in files_list:
        with open(file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("ID"):
                    print("ID={0}".format(counter))
                    counter = counter + 1
                else:
                    print(line)


def main():
    if len(sys.argv) < 2:
        print("Usage: merge_msalign.py msalign_file_1[,msalign_file_2...]\n\n"
              "Merges multiple .msalign spectra files into one",
              file=sys.stderr)
        return 1

    merge(sys.argv[1:])
    return 0


if __name__ == "__main__":
    main()
