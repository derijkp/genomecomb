#!/bin/sh
#$ -S /bin/bash
#$ -V
#$ -cwd
cd $1
shift
exec "$@"
