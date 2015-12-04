#!/usr/bin/env python

#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

from __future__ import print_function
import shutil
import sys
import os
import argparse
import subprocess
from threading import Thread
from time import sleep

spectrogene_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, spectrogene_root)

LIB_DIR = "lib"
lib_absolute = os.path.join(spectrogene_root, LIB_DIR)
sys.path.insert(0, lib_absolute)
sys.path.insert(0, spectrogene_root)
os.environ["PATH"] += os.pathsep + lib_absolute


from spectrogene.process_proteome import ProteomeProcessor
from spectrogene.main import (_run_toppic, _read_spectra, _merge_toppic_tables,
                              _write_spectra, _split_strings_list, _copy_html)


def _run_parallel_proteome(input_proteome, input_spectra, work_dir, num_proc):
    """
    Runs toppic in parallel threads
    """
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    work_dir = os.path.abspath(work_dir)

    input_spectra = os.path.abspath(input_spectra)
    input_proteome = os.path.abspath(input_proteome)

    spectra_text = _read_spectra(input_spectra)
    spec_splitted = _split_strings_list(spectra_text, num_proc)
    threads = []
    output_files = []

    for i in range(num_proc):
        inst_name = "part_{0}".format( i)
        inst_workdir = os.path.join(work_dir, inst_name)
        if os.path.isdir(inst_workdir):
            shutil.rmtree(inst_workdir)
        os.mkdir(inst_workdir)

        inst_spec = os.path.join(inst_workdir, "spectra.msalign")
        _write_spectra(spec_splitted[i], inst_spec)

        inst_prot = os.path.join(inst_workdir, "proteome.fasta")
        shutil.copy2(input_proteome, inst_prot)

        print("Running {0} TopPic instance".format(inst_name))
        thread = Thread(target=_run_toppic,
                        args=(inst_prot, inst_spec, inst_workdir))
        thread.start()
        threads.append(thread)

        out_file = os.path.join(inst_workdir, "spectra.OUTPUT_TABLE")
        output_files.append(out_file)

    for t in threads:
        t.join()

    return output_files


def main():
    parser = argparse.ArgumentParser(description="Runs SpectroGene in proteome"
                                     " mode (for benchmarking purposes)",
                                     formatter_class= \
                                        argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("prot_file", help="path to proteome file in "
                        "FASTA fromat")
    parser.add_argument("prot_coords", help="a file with protein "
                                                     "coordinates")
    parser.add_argument("genome_file", help="path to genome in FASTA fromat")
    parser.add_argument("spectra_file", help="path to spectra file in "
                        "MsAlign format")
    parser.add_argument("-o", "--output_dir", help="output directory",
                        dest="output_dir", default="specttrogene-out")
    parser.add_argument("-p", "--num_proc", dest="num_proc", type=int,
                                 help="number of processes", default="1")
    parser.add_argument("-e", "--evalue", dest="e_value", type=float,
                        help="E-value threshold for PrSM", default="0.01")
    args = parser.parse_args(sys.argv[1:])


    merged_output = os.path.join(args.output_dir, "toppic_results.txt")
    if not os.path.isfile(merged_output):
        out_files = _run_parallel_proteome(args.prot_file, args.spectra_file,
                                          args.output_dir, args.num_proc)
        _merge_toppic_tables(out_files, open(merged_output, "w"))
    else:
        print("Using TopPic results from the previous run")

    proc = ProteomeProcessor(args.e_value, args.genome_file,
                             args.prot_coords)
    out_prsms = os.path.join(args.output_dir, "prsms.txt")
    proc.process(merged_output)

    out_prsms = os.path.join(args.output_dir, "prsms.txt")
    proc.output_prsms(out_prsms)

    out_orfs = os.path.join(args.output_dir, "orf_clusters.txt")
    proc.print_orfs(out_orfs)

    html_dir = os.path.join(args.output_dir, "prsm_html")
    _copy_html(proc.prsms, html_dir)


if __name__ == "__main__":
    main()
