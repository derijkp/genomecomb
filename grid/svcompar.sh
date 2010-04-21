#!/bin/sh
#$ -S /bin/bash
#$ -V
#$ -cwd
cg svcompare GS103/GS103-$1-paired-sv.tsv GS102/GS102-$1-paired-sv.tsv
