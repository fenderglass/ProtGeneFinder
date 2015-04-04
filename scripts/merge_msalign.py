#!/usr/bin/env python

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
    merge(sys.argv[1:])


if __name__ == "__main__":
    main()
