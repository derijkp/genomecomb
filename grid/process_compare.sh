#!/bin/sh
#$ -S /bin/bash
#$ -V
#$ -cwd
export dbdir=/mnt/extra/CompleteGenomics/refseq
cg process_compare $dbdir $1
