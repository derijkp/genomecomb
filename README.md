GenomeComb
========== 
A program for the analysis of genome and transciptome sequencing data
Copyright VIB and University of Antwerp

Purpose
-------
Genomecomb is a package designed to analyze, combine, annotate
and query genome as well as transcriptome sequencing data.

It provides pipelines to analyse data sets made with different sequencing
technologies (Illumina, nanopore) using a (choosable) set of aligners,
variant callers, structural variant callers, isoform callers, gene
counters, etc. and to extensively annotate the results.

It also provides a very flexible tool to query or summarise the results, and
includes a graphical table browser that can handle browsing (and querying)
multigigabyte tables with millions of rows

Other included tools can be used for region operations (including making a
queryable multisample regionfile), validation (primer design and
checking), genetics analysis of sequencing data (homwes, plink export),
conversion (between different formats, extended liftover) and reporting
(depth histograms, pipeline results, gender prediction)

Availability
------------
Genomecomb is available on github ([https://github.com/derijkp/genomecomb](https://github.com/derijkp/genomecomb))
General information and extensive documentation can be found on the pages website
[https://derijkp.github.io/genomecomb](https://derijkp.github.io/genomecomb)

System requirements and Installation
------------------------------------
Binary packages are available for Linux. Genomecomb itself is distributed as a
portable application directory: A self-contained directory with the
genomecomb executable (cg) and all needed depencies compiled in a way they
should work on all (except very ancient) Linux systems.

Installation of the package is as simple as downloading the 
[archive](https://github.com/derijkp/genomecomb/releases/download/0.109.0/genomecomb-0.109.0-Linux-x86_64.tar.gz)
and unpacking it somewhere. You can either call the executable (cg) directly
from the directory, or put a soft-link to it somehwere in the path. (The
executable itself must stay in the application directory to work.)

Most analyses need a reference genome and accompanying annotation
databases. These can be downloaded separately for a number of species as
described in the [genomecomb installation documentation](https://derijkp.github.io/genomecomb/install.html)

Some of the externally developed software used by genomecomb are included
in the distribution (e.g. samtools). However most external software (such
as e.g specific variant callers) needs to be installed on the system.
Portable application directories of these programs (in versions that were
tested with genomecomb) are made available on the the [genomecomb
website](https://derijkp.github.io/genomecomb/install.html). These can be
made generally available (put in the PATH) or to genomecomb only by
putting them in the directory extra in the genomecomb app directory.

Installation quickstart example:
```
cd ~/bin
wget https://genomecomb.bioinf.be/download/genomecomb-0.109.0-Linux-x86_64.tar.gz
tar xvzf genomecomb-0.109.0-Linux-x86_64.tar.gz

# optional if ~/bin is not in your PATH already
export PATH=$HOME/bin:$PATH
# or make a softlink to something in the PATH already
ln -s ~/bin/genomecomb-0.109.0-Linux-x86_64/cg /usr/local/bin/cg

# Install the hg38 reference databases and the software needed to run the ont and srs (short read) presets.
cg install hg38 hg38-cadd hg38-minimap2 srs ont

```

Documentation
-------------
A general overview of how to use the program is given in the help.
You can get help from the program itself using the command
```
cg help
```
You can also find extensive documentation (including all help pages) on 
the [genomecomb documentation Pages](https://derijkp.github.io/genomecomb/intro.html)

Bugs
----
If you are having problems with the program contact me. I will
do my best to get it fixed. Please report any bugs you have
found. If possible, state your machine's hardware and software
configurations. Sending me a full description of the
circumstances in which the bug occurs, possibly with the data it
happened on, will help me tracking down a bug. If you have any
suggestions, you can also make them to me.

License
-------
The use of this application is governed by the GPL (license.txt).

We have used a very early version of GenomeComb in the analysis of whole
genome sequences of monozygotic twin genomes, tumor-normal genomes and
publicly available HapMap genomes. A citable publication about this study
has appeared in Nature Biotch:

* Reumers, J*, De Rijk, P*, Zhao, H, Liekens, A, Smeets, D, Cleary, J, Van Loo, P, Van Den Bossche, M, Catthoor, K, Sabbe, B, Despierre, E, Vergote, I, Hilbush, B, Lambrechts, D and Del-Favero, J ;
Optimized filtering reduces the error rate in detecting genomic variants by short-read sequencing.
Nature biotechnology, 30, 61-88
pubmed 22178994

For the homwes tool available in genomecomb the folowing paper should be referenced:

* Kancheva,D, Atkinson,D, De Rijk,P, Zimon,M, Chamova,T, Mitev,V, Yaramis,A, Maria Fabrizi,G, Topaloglu,H, Tournev,I, Parma, Y, Battaloglu, E, Estrada-Cuzcano, A, Jordanova, A (2015)
Novel mutations in genes causing hereditary spastic paraplegia and Charcot-Marie-Tooth neuropathy identified by an optimized protocol for homozygosity mapping based on whole-exome sequencing.
Genet. Med., 10.1038

How to contact me
-----------------

Peter De Rijk
VIB - UAntwerp Center for Molecular Neurology, Neuromics Support Facility - Bioinformatics
University of Antwerp
Universiteitsplein 1
B-2610 Antwerpen, Belgium

tel.: +32-03-265.10.40
E-mail: Peter.DeRijk@uantwerpen.be
