#!/bin/sh
#$ -S /bin/bash
#$ -V
#$ -cwd
export multicompar=/complgen/multicompar/compar.tsv
cg compare_mpvt $multicompar $1
