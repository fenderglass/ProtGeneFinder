SpectroGene
===========

Version 1.0

Description
-----------

SpectroGene is a tool for top-down protein identification
using unannotated bacterial genome. SpectroGene could be used for disocvery
of previously unannotated genes as well as correction of known genes coordinates
(especially, Start codons). Almost all PrSMs which are found with a
standard search with a known proteome could be also identified with SpectroGene
without proteome reference.


Availability
------------

Currently only Linux x86-64 systems are supported (because of TopPIC dependency).


Requirements
------------

* Python 2.7
* TopPIC (included) [http://proteomics.informatics.iupui.edu/software/toppic/]
* Biopython [biopython.org]


Install
-------
SpectroGene is written on Python and requires no installation. 
However, it depends on some extra Python modules listed in "Requirements" section.
These dependencies could be installed via your system's package manager or "pip"
(e.g. "pip install biopython").


Quick Usage
-----------
    usage: spectrogene.py [-h] [-o OUTPUT_DIR] [-p NUM_PROC] [-e E_VALUE]
                          [--version]
                          genome_file spectra_file

    SpectroGene v1.0: Genome annotation using top-down mass spectra

    positional arguments:
      genome_file           path to genome file in FASTA fromat
      spectra_file          path to spectra file in MsAlign format

    optional arguments:
       -h, --help            show this help message and exit
       -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                            output directory (default: spectrogene-out)
      -p NUM_PROC, --num_proc NUM_PROC
                            number of processes (default: 1)
      -e E_VALUE, --evalue E_VALUE
                            E-value threshold for PrSM (default: 0.01)
      --version             show program's version number and exit


Detailed Instructions
---------------------



Output formats
--------------



Authors
-------
* Mikhail Kolmogorov (UCSD)


Citation
--------
* Mikhail Kolmogorov, Xiaowen Liu and Pavel Pevzner. "SpectroGene: a tool for proteogenomic annotations using top-down spectra"


Contacts
--------
Please report any bugs directly to the issue tracker of this project.
Also, you can send your feedback at fenderglass@gmail.com