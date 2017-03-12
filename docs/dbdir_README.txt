Reference and annotation databases for genomecomb
=================================================
version: 0.98.0

This directory contains reference and annotation databases for genomecomb. 
It is meant to used as the "refdb" directory in many of the genomecomb commands.

Reference genome
----------------
The genome is present as an ifas file, named genome_build.ifas where build is
replaced by the reference sequence build id (e.g. hg19).
The ifas format is completely compatible with the fasta format, 
but adds a few restrictions
* Each sequence is on one continuous line (no returns allowed).
* The sequences are sorted according to chromosome name (using natural sort).

The file genome_build.ifas.index contains a list of the chromosome names
and the start position of their sequence (in the ifas file).

Further genome indexes can be present, e.g.
* genome_build.ssa: sussinct suffix array index of the genome used by several genomecomb tools (e.g. primercheck)
* genome_build.ifas.bwa: directory containing the files necesary for using the bwa aligner

Annotation databases
--------------------
All annotation files/databases in the formats accepted by "cg annotate" in
this directory will be used to annotate variants, etc. when used as
"refdb" dir. Each annotation databases has a filename according to the
following structure: type_build_name.ext
The extension indicates the type of file (tsv or bcol).
The part of the filename before the first _ indicates the type of annotation:
* reg: region file that will be used to annotate regions that it overlaps with
* var: a variations file to annotate variants that must match to the given location, type and allele
* gene: gene files (in genepred-like format) are used to annotate variants with
the effects they have on the genes.
* mir: a miRNA gene file
* bcol: bcol databases are used to annotate positions (e.g. snps) with a given value. Database 
files are in the [[bcol]] format (also extension bcol).

Together with the annotation file is often an info file (same name but
with extension info) containig information about the database and an opt
file (same name but with extension opt) that contains settings for the
annotation

directory extra
---------------
Annotation files in this subdirectory extra are not used to annotate files
by default because they are only useful as annotation in specific cases or
mainly useful for analysis other than annotation. e.g.
extra/reg_build_fullgenome.tsv (a region file of the full genome) is not a
useful annotation, but can be useful to subtract other region files from
to see what is left from the genome.
