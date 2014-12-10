#!/usr/bin/env python

from __future__ import print_function
import shutil
import os
import sys
import subprocess
from threading import Thread
from time import sleep
import argparse

from make_proteome import make_proteome

MSALIGN_DIR = os.environ["MSALIGN_DIR"]
MSALIGN_CMD = ["java", "-Xmx12G", "-classpath", "jar/*:",
               "edu.ucsd.msalign.align.console.MsAlignPipeline"]

WINDOWS = [500, 230, 110, 60]

def run_instance(proteome_file, spectrum_file, work_dir):
    assert os.path.isdir(work_dir)
    os.mkdir(os.path.join(work_dir, "msinput"))
    os.mkdir(os.path.join(work_dir, "msoutput"))
    os.mkdir(os.path.join(work_dir, "html"))
    os.mkdir(os.path.join(work_dir, "xml"))
    shutil.copytree(os.path.join(MSALIGN_DIR, "etc"),
                    os.path.join(work_dir, "etc"))
    shutil.copytree(os.path.join(MSALIGN_DIR, "xsl"),
                    os.path.join(work_dir, "xsl"))

    input_config = os.path.join(MSALIGN_DIR, "msinput", "input.properties")
    shutil.copy2(proteome_file, os.path.join(work_dir, "msinput", "prot.fasta"))
    shutil.copy2(spectrum_file, os.path.join(work_dir, "msinput",
                                             "spectra.msalign"))
    shutil.copy2(input_config, os.path.join(work_dir, "msinput"))

    subprocess.check_call(MSALIGN_CMD + [work_dir],
                          stdout=open(os.devnull, "w"))


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
    os.chdir(MSALIGN_DIR)

    spectras_text = read_spectrum_file(input_spectrum)
    spec_splitted = split_n(spectras_text, num_proc)
    threads = []

    for i in range(num_proc):
        inst_name = "part_{0}".format( i)
        inst_workdir = os.path.join(work_dir, inst_name)
        if os.path.isdir(inst_workdir):
            shutil.rmtree(inst_workdir)
        os.mkdir(inst_workdir)

        inst_spec = os.path.join(inst_workdir, "spectra.msalign")
        write_spectras(spec_splitted[i], inst_spec)

        print("Running {0} instance".format(inst_name))
        thread = Thread(target=run_instance,
                        args=(input_proteome, inst_spec, inst_workdir))
        thread.start()
        threads.append(thread)

    for t in threads:
        t.join()


def run_parallel_genome(input_genome, input_spectrum, work_dir, num_proc):
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    work_dir = os.path.abspath(work_dir)

    input_spectrum = os.path.abspath(input_spectrum)
    input_genome = os.path.abspath(input_genome)
    os.chdir(MSALIGN_DIR)

    for window in WINDOWS:
        print("Running with window size", window)
        window_dir = os.path.join(work_dir, "genome_run_" + str(window))
        if not os.path.isdir(window_dir):
            os.mkdir(window_dir)

        print("Creating proteome")
        prot_file_1, prot_file_2 = make_proteome(input_genome, window,
                                                 window_dir)
        print("Reading spectrum")
        spectras_text = read_spectrum_file(input_spectrum)

        proc_per_half = int(num_proc / 2)
        spec_splitted = split_n(spectras_text, proc_per_half)
        threads = []

        def run_for_half(inst_pref, inst_prot):
            for i in range(proc_per_half):
                inst_name = "{0}_{1}".format(inst_pref, i)
                inst_workdir = os.path.join(window_dir, inst_name)
                if os.path.isdir(inst_workdir):
                    shutil.rmtree(inst_workdir)
                os.mkdir(inst_workdir)

                inst_spec = os.path.join(inst_workdir, "spectra.msalign")
                write_spectras(spec_splitted[i], inst_spec)

                print("Running {0} instance".format(inst_name))
                thread = Thread(target=run_instance,
                                args=(inst_prot, inst_spec, inst_workdir))
                thread.start()
                threads.append(thread)

        run_for_half("noshift", prot_file_1)
        run_for_half("halfshift", prot_file_2)

        for t in threads:
            t.join()


def main():
    parser = argparse.ArgumentParser(description="Run MSAlign on genome/proteome")
    subparsers = parser.add_subparsers(help="sub-command help", dest="mode")

    parser_proteome = subparsers.add_parser("proteome", help="proteome mode")
    parser_proteome.add_argument("prot_file", help="path to proteome file")
    parser_proteome.add_argument("spectrum_file", help="path to spectrum file")
    parser_proteome.add_argument("output_dir", help="output directory")
    parser_proteome.add_argument("num_proc", help="number of threads")

    parser_genome = subparsers.add_parser("genome", help="genome mode")
    parser_genome.add_argument("genome_file", help="path to genome file")
    parser_genome.add_argument("spectrum_file", help="path to spectrum file")
    parser_genome.add_argument("output_dir", help="output directory")
    parser_genome.add_argument("num_proc", help="number of threads")

    args = parser.parse_args(sys.argv[1:])

    if args.mode == "genome":
        run_parallel_genome(args.genome_file, args.spectrum_file,
                            args.output_dir, int(args.num_proc))
    else:
        run_parallel_proteome(args.prot_file, args.spectrum_file,
                              args.output_dir, int(args.num_proc))

    return 0


if __name__ == "__main__":
    main()
