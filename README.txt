GenomeComb
========== A program for the analysis of complete genome data
           Copyright 2010 to 2021 VIB and Universiteit Antwerpen

Purpose
-------
Genomecomb is a package designed to analyze, combine, annotate
and query whole genome, exome or targetted sequencing data.

It provides pipelines to analyse data sets made with different sequencing
technologies (Illumina, nanopore) using a (choosable) set of aligners,
variant callers, structural variant callers, etc. and to extensively
annotate the results.

It also provides a very flexible tool to query or summarise the results, and
includes a graphical table browser that can handle browsing (and querying)
multigigabyte tables with millions of rows

Other included tools can be used for region operations (including making a
queryable multisample regionfile), validation (primer design and
checking), genetics analysis of sequencing data (homwes, plink export),
conversion (between different formats, extended liftover) and reporting
(depth histograms, pipeline results, gender prediction)

System requirements and Installation
------------------------------------
Updates and more information can be found on the website
http://genomecomb.sourceforge.net/

Binary packages are available for Linux. The application 
uses the concept of an application directory: Everything needed to run
the application is in one directory: installation is as simple as moving
the directory to the place you want it (unpacking it, if it is compressed).
Uninstallation is done by removing the directory.
The executable (cg) must stay in the application directory to work.
If you want to be able to start it from another location (desktop, bin),
make a link to it.

A general overview of how to use the program is given in the help.
You can get help from the program itself using the command
cg help

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
The use of this application is governed by the accompanying 
license (license.txt).
If you have used the program to obtain results, you should cite 
the following paper:

Reumers, J*, De Rijk, P*, Zhao, H, Liekens, A, Smeets, D, Cleary, J, Van Loo, P, Van Den Bossche, M, Catthoor, K, Sabbe, B, Despierre, E, Vergote, I, Hilbush, B, Lambrechts, D and Del-Favero, J
Optimized filtering reduces the error rate in detecting genomic variants by short-read sequencing.
Nature biotechnology, 30, 61-88
pubmed 22178994

How to contact me
-----------------

Peter De Rijk
VIB - UAntwerp Center for Molecular Neurology, Neuromics Support Facility - Bioinformatics
University of Antwerp
Universiteitsplein 1
B-2610 Antwerpen, Belgium

tel.: +32-03-265.10.40
E-mail: Peter.DeRijk@uantwerpen.vib.be
