portable hap.py
===============
made using script make_extra_happy.sh from genomecomb (Peter de Rijk, Peter.DeRijk@uantwerpen.vib.be)

This is a portable directory version of hap.py (https://github.com/Illumina/hap.py).
It comes with all dependencies needed in one directory (compiled for wide compatibility), so
you can run it on about any linux distribution (except very old ones) without installation:
just start the hap.py "executable" in the directory.

The executable must stay in the portable directory to work.
If you want the command to be avialable without specifying the path, you have to include
the portable directory in the PATH, or make a soft link to happy-0.3.14/hap.py in a (bin) directory
that is in the PATH.

This version was created using conda using the holy build box
then extracted in a directory using conda pack.
The hap.py "executable" placed in the dir is a short shell scripts that adapts the environment
so the portable dir resources are used and then runs hap.py