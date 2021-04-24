#!/bin/bash
# This script will run hap.py (with the given parameters) in its own environment
# (with matching python, libs, etc.)
# It will find the real (environment) directory hap.py is in by following softlinks
# It will change the PATH to include the <environment dir>/bin
# It will change the LD_LIBRARY_PATH to include the <environment dir>/lib
# and will then run hap.py with the given parameters
script="$(readlink -f "$0")"
dir="$(dirname "$script")"
PATH=$dir/bin:$PATH
LD_LIBRARY_PATH=$dir/lib:$LD_LIBRARY_PATH
hap.py ${1+"$@"}
