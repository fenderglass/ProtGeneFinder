#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

"""
The main SpectroGene module. It defines top-level logic of the program
"""

from __future__ import print_function
import shutil
import os
import sys
import subprocess
from threading import Thread
from time import sleep
import argparse

from spectrogene.make_proteome import make_proteome
from spectrogene.process_genome import GenomeProcessor
from spectrogene.process_proteome import ProteomeProcessor
from spectrogene.__version__ import __version__

TOPPIC_BIN = "toppic"
WINDOWS = [500, 200, 50]


def _run_toppic(proteome_file, spectra_file, work_dir):
    """
    Invokes TopPic binary
    """
    cmdline = [TOPPIC_BIN, proteome_file, spectra_file, "--max-ptm", "500",
               "--generating-function", "--cutoff-value", "100"]
    subprocess.check_call(cmdline, stdout=open(os.devnull, "w"))


def _merge_toppic_tables(tables, out_stream):
    """
    Merges multiple TopPic output tables into a single one
    """
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


def _read_spectra(filename):
    """
    Reads spectra file in msalign format
    """
    spectra = []
    curent_spec = ""
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                spectra.append(curent_spec)
                curent_spec = ""
                continue

            curent_spec += line + "\n"

    if curent_spec:
        spectra.append(curent_spec)

    return spectra


def _write_spectra(spectra_strings, out_file):
    """
    Writes a list of spectra texts into file
    """
    with open(out_file, "w") as f:
        for ss in spectra_strings:
            f.write(ss + "\n")


def _split_strings_list(lst, n):
    """
    Splits a list of strings into n chunks of roughly equal size
    """
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


def _run_parallel_genome(input_genome, input_spectra, work_dir, num_proc):
    """
    Generates ORFeome and runs toppic in parallel threads
    given genome and spectra files
    """
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    work_dir = os.path.abspath(work_dir)

    input_spectra = os.path.abspath(input_spectra)
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
        print("Reading spectra")
        spectra_text = _read_spectra(input_spectra)
        spec_splitted = _split_strings_list(spectra_text, num_proc)
        threads = []

        for i in range(num_proc):
            inst_name = "part_" + str(i)
            inst_workdir = os.path.join(window_dir, inst_name)
            if os.path.isdir(inst_workdir):
                shutil.rmtree(inst_workdir)
            os.mkdir(inst_workdir)

            inst_spec = os.path.join(inst_workdir, "spectra.msalign")
            _write_spectra(spec_splitted[i], inst_spec)

            inst_prot = os.path.join(inst_workdir, "proteome.fasta")
            shutil.copy2(prot_file, inst_prot)

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


def _copy_html(prsms, html_dir):
    """
    Copies html annotations generated by TopPic
    """
    #TODO: copy necessary resource files
    if not os.path.isdir(html_dir):
        os.mkdir(html_dir)

    for prsm in prsms:
        html_name = os.path.join(html_dir, "spec{0}.html".format(prsm.spec_id))
        shutil.copy2(prsm.html, html_name)


def main():
    parser = argparse.ArgumentParser(description="SpectroGene v{0}: Genome "
                                    "annotation using top-down mass spectra"
                                     .format(__version__),
                                     formatter_class= \
                                        argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("genome_file",
                        help="path to genome file in FASTA fromat")
    parser.add_argument("spectra_file",
                        help="path to spectra file in MsAlign format")
    parser.add_argument("-o", "--output_dir", dest="output_dir",
                        default="spectrogene-out", help="output directory")
    parser.add_argument("-p", "--num_proc", dest="num_proc", type=int,
                               help="number of processes", default="1")
    parser.add_argument("-e", "--evalue", dest="e_value", type=float,
                            help="E-value threshold for PrSM",
                            default="0.01")
    parser.add_argument("--version", action="version", version=__version__)

    args = parser.parse_args(sys.argv[1:])

    merged_output = os.path.join(args.output_dir, "toppic_results.txt")
    if not os.path.isfile(merged_output):
        out_files = _run_parallel_genome(args.genome_file, args.spectra_file,
                                        args.output_dir, args.num_proc)
        _merge_toppic_tables(out_files, open(merged_output, "w"))
    else:
        print("Using TopPic results from the previous run")

    proc = GenomeProcessor(args.e_value, args.genome_file)
    proc.process(merged_output)

    out_prsms = os.path.join(args.output_dir, "prsms.txt")
    proc.output_prsms(out_prsms)

    out_orfs = os.path.join(args.output_dir, "orf_clusters.txt")
    proc.print_orfs(out_orfs)

    html_dir = os.path.join(args.output_dir, "prsms_html")
    _copy_html(proc.prsms, html_dir)

    return 0
