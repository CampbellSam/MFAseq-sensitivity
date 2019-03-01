# README #

This repository is to share code used in the analysis of _Leishmania_ origins of DNA replication in the paper.

>Genome-wide mapping reveals single-origin chromosome replication in _Leishmania_, a eukaryotic microbe
>Paper authors: Catarina A. Marques, Nicholas J. Dickens, Daniel Paape, Samantha J. Campbell and Richard McCulloch

Software authors are detailed in contributors.txt
### What is in this repository? ###

* generate_simulated_origins.py (version 1.0) - this script was used to simulate origins of DNA replication at strand switch regions (SSR) using existing data and a known origin.
* readBed.py - this script was used to measure the ratios across whole chromosomes, SSR and in a window around the putative origin
* measure_peak_heights.sh - an example of batch usage of the script, used on a batch of bed files with ratios (eS/G2, etc)
* LICENSE.txt is a copy of the Gnu Public License version 3, all software and documents in this repository are distributed under this license whether explicity described inside it or not.

### Copyright ###
All software in this repository is copyright 2015, Nicholas J. Dickens and Samantha J. Campbell, address Wellcome Trust Centre for Molecular Parasitology, University of Glasgow, 120 University Place, Glasgow, G12 8TA, UK.

### Dependencies ###

There are very few dependencies and these scripts require no special installation.

* Python >=2.7
* pysam (0.8.3 recommended)
* numpy >= 1.8.2
* pybedtools (required by readBed.py)

* bedtools installed in your path (required by pybedtools)

### Who do I talk to? ###

For further information about the scripts in this repository please contact Nick Dickens by [e-mail](mailto:nick.dickens@glasgow.ac.uk) or [twitter](http://twitter.com/WTCMPbix)

