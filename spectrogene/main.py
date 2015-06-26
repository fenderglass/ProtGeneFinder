from __future__ import print_function
import shutil
import os
import sys
import subprocess
from threading import Thread
from time import sleep
import argparse

from spectrogene.make_proteome import make_proteome
from spectrogene.merge_tables import merge_tables
from spectrogene.process_genome import GenomeProcessor
from spectrogene.process_proteome import ProteomeProcessor

TOPPIC_BIN = "toppic"
WINDOWS = [500, 200, 50]

def run_instance(proteome_file, spectrum_file, work_dir):
    cmdline = [TOPPIC_BIN, proteome_file, spectrum_file, "--max-ptm", "500",
               "--generating-function", "--cutoff-value", "100"]
    subprocess.check_call(cmdline, stdout=open(os.devnull, "w"))


def read_spectrum_file(filename):
    spectras = []
    curent_spec = ""
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                spectras.append(curent_spec)
                curent_spec = ""
                continue

            curent_spec += line + "\n"

    if curent_spec:
        spectras.append(curent_spec)

    return spectras


def write_spectras(spectra_strings, out_file):
    with open(out_file, "w") as f:
        for ss in spectra_strings:
            f.write(ss + "\n")


def split_n(lst, n):
    total_len = sum(map(len, lst))
    avg = total_len / float(n)
    out_lst = []
    last = 0.0
    lst_cur = 0

    while len(out_lst) < n:
        slice = []
        slice_len = 0.0

        while (last < avg * (len(out_lst) + 1) and
               lst_cur < len(lst)):
            slice.append(lst[lst_cur])
            last += len(slice[-1])
            lst_cur += 1

        out_lst.append(slice)

    assert len(out_lst) == n
    return out_lst


def run_parallel_proteome(input_proteome, input_spectrum, work_dir, num_proc):
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    work_dir = os.path.abspath(work_dir)

    input_spectrum = os.path.abspath(input_spectrum)
    input_proteome = os.path.abspath(input_proteome)

    spectras_text = read_spectrum_file(input_spectrum)
    spec_splitted = split_n(spectras_text, num_proc)
    threads = []
    output_files = []

    for i in range(num_proc):
        inst_name = "part_{0}".format( i)
        inst_workdir = os.path.join(work_dir, inst_name)
        if os.path.isdir(inst_workdir):
            shutil.rmtree(inst_workdir)
        os.mkdir(inst_workdir)

        inst_spec = os.path.join(inst_workdir, "spectra.msalign")
        write_spectras(spec_splitted[i], inst_spec)

        inst_prot = os.path.join(inst_workdir, "proteome.fasta")
        shutil.copy2(input_proteome, inst_prot)

        print("Running {0} instance".format(inst_name))
        thread = Thread(target=run_instance,
                        args=(inst_prot, inst_spec, inst_workdir))
        thread.start()
        threads.append(thread)

        out_file = os.path.join(inst_workdir, "spectra.OUTPUT_TABLE")
        output_files.append(out_file)

    for t in threads:
        t.join()

    return output_files


def run_parallel_genome(input_genome, input_spectrum, work_dir, num_proc):
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    work_dir = os.path.abspath(work_dir)

    input_spectrum = os.path.abspath(input_spectrum)
    input_genome = os.path.abspath(input_genome)
    output_files = []

    for window in WINDOWS:
        print("Running with window size", window)
        window_dir = os.path.join(work_dir, "genome_run_" + str(window))
        if not os.path.isdir(window_dir):
            os.mkdir(window_dir)

        print("Creating proteome")
        prot_file = os.path.join(window_dir, "proteome.fasta")
        make_proteome(input_genome, window, prot_file)
        print("Reading spectrum")
        spectras_text = read_spectrum_file(input_spectrum)
        spec_splitted = split_n(spectras_text, num_proc)
        threads = []

        for i in range(num_proc):
            inst_name = "part_" + str(i)
            inst_workdir = os.path.join(window_dir, inst_name)
            if os.path.isdir(inst_workdir):
                shutil.rmtree(inst_workdir)
            os.mkdir(inst_workdir)

            inst_spec = os.path.join(inst_workdir, "spectra.msalign")
            write_spectras(spec_splitted[i], inst_spec)

            inst_prot = os.path.join(inst_workdir, "proteome.fasta")
            shutil.copy2(prot_file, inst_prot)

            print("Running {0} instance".format(inst_name))
            thread = Thread(target=run_instance,
                            args=(inst_prot, inst_spec, inst_workdir))
            thread.start()
            threads.append(thread)

            out_file = os.path.join(inst_workdir, "spectra.OUTPUT_TABLE")
            output_files.append(out_file)

        for t in threads:
            t.join()

    return output_files


def main():
    parser = argparse.ArgumentParser(description="Genome annotation using "
                                               "top-down mass spectrometry")

    subparsers = parser.add_subparsers(help="sub-command help", dest="mode")

    parser_proteome = subparsers.add_parser("proteome", help="proteome mode")
    parser_proteome.add_argument("prot_file", help="path to proteome file")
    parser_proteome.add_argument("prot_coords", help="a file with protein "
                                                     "coordinates")
    parser_proteome.add_argument("genome_file", help="path to fasta file")
    parser_proteome.add_argument("spectrum_file", help="path to spectrum file")
    parser_proteome.add_argument("output_dir", help="output directory")
    parser_proteome.add_argument("-p", "--num_proc", dest="num_proc", type=int,
                                 help="number of processes", default="1")
    parser_proteome.add_argument("-e", "--evalue", dest="e_value", type=float,
                            help="E-value threshold for PrSM [default 0.01]",
                            default="0.01")

    parser_genome = subparsers.add_parser("genome", help="genome mode")
    parser_genome.add_argument("genome_file", help="path to genome file")
    parser_genome.add_argument("spectrum_file", help="path to spectrum file")
    parser_genome.add_argument("output_dir", help="output directory")
    parser_genome.add_argument("-p", "--num_proc", dest="num_proc", type=int,
                               help="number of processes", default="1")
    parser_genome.add_argument("-e", "--evalue", dest="e_value", type=float,
                            help="E-value threshold for PrSM [default 0.01]",
                            default="0.01")

    args = parser.parse_args(sys.argv[1:])

    if args.mode == "genome":
        #out_files = run_parallel_genome(args.genome_file, args.spectrum_file,
        #                                args.output_dir, args.num_proc)
        proc = GenomeProcessor(args.e_value, args.genome_file)
    else:
        #out_files = run_parallel_proteome(args.prot_file, args.spectrum_file,
        #                                  args.output_dir, args.num_proc)
        proc = ProteomeProcessor(args.e_value, args.genome_file,
                                 args.prot_coords)

    merged_output = os.path.join(args.output_dir, "toppic_merged.txt")
    #merge_tables(out_files, open(merged_output, "w"))
    out_prsms = os.path.join(args.output_dir, "prsms.sg")
    out_orfs = os.path.join(args.output_dir, "orf_clusters.sg")
    proc.process_and_output(merged_output, out_prsms)
    proc.print_orfs(out_orfs)
    #proc.copy_html(args.output_dir)

    return 0
