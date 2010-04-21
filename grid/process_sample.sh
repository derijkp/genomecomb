#!/bin/sh
#$ -S /bin/bash
#$ -V
#$ -cwd
export dbdir=/mnt/extra/CompleteGenomics/refseq
cg process_sample $1 $dbdir
