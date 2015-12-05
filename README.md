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
* Biopython [http://biopython.org]


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


Running Instructions
--------------------

SpectroGene takes as input:

* Genome file in FASTA format
* Deconvoluted top-down spectra in msalign format (similar to MGF)

A typical pipeline for raw spectra (such as from Thermo Scientific™ Orbitrap™ 
mass spectrometers) would include centroiding and deconvolution.

First, raw spectra should be centroided and converted to mzXML format,
which could be done using msconvert utility 
[available for Windows only: http://proteowizard.sourceforge.net/tools.shtml].
Here is an example command line:

    msconvert.exe --filter "peakPicking true 1-" --mzXML --zlib spectra.raw

Then, the resulting mzXML file should be deconvoluted using MsDeconv
software [http://bix.ucsd.edu/projects/msdeconv/]:

    java -jar msdeconv.jar spectra.mzXML -o msalign

Finally, you can run SpectroGene on resulting .msalign file:

    ./spectrogene.py genome.fasta spectra.msalign


Output
------

SpectroGene provides three kinds of output:

* PrSMs mapped on the genome sequence
* HTML visualisation of PrSMs (done by TopPIC)
* ORF clusters of identified proteoforms

PrSMs are output into a file called "prsms.txt" in the output directory.
It is formatted as CSV table and have self-explanatory columns.
A graphical representations of all identified spectra in HTML
fromat could be found in "html_prsms" directory.

Also, raw TopPIC output with additional spectra information is available 
in "toppic_merged.txt" file.

The information about identified ORF cluster is output into "orf_clusters.txt"
file. For each cluster multiple statistics are provided (such as positions,
number of spectra/proteoforms etc.). Below the statistics, proteoforms aligned
on the ORF sequence are shown. Putative Start codons are marked with "v"
symbol on the top of the ORF sequence. Each proteoform might have multiple corresponding
spectra. In that case, we choose the spectra with the lowest P-value as representative.
Representative spectra ids and P-values are shown on the right side of the proteoforms.
Also, post-translational modifications are shown above the corresponding
proteoforms (with numbers corresponding to a modification weight).


Useful scripts
--------------

SpectrGene distribution also contain some additional scripts
which might be useful (located in "scripts" directory):

**merge_msalign.py**

Merges multiple msalign spectra files into one.

**statistics.py**

Provides some extra statistics on the SpectroGene run.

**compare.py**

Compares two SpectrGene runs (for benchmarking purposes)

**proteome_run.py**

Runs SpectroGene with a known proteome (also for benchmarking, see the paper for details).


License
-------
The program is distributed under BSD license. See LICENSE file for details.


Citation
--------
* Mikhail Kolmogorov, Xiaowen Liu and Pavel Pevzner. "SpectroGene: a tool for proteogenomic annotations using top-down spectra"


Authors
-------
* Mikhail Kolmogorov (UCSD)


Contacts
--------
Please report any bugs directly to the issue tracker of this project.
Also, you can send your feedback at fenderglass@gmail.com