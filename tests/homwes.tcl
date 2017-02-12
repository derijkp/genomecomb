#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]

test homwes {basic} {
	file delete -force $::bigtestdir/tmp/homwes
	file mkdir $::bigtestdir/tmp/homwes
	set varfile $::bigtestdir/ori/exomes_yri/compar/annot_compar-yri_exome_test.tsv.lz4
	set dbdir $::bigtestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile NA19240 $::bigtestdir/tmp/homwes/homwes-out.tsv
	checkdiff -qr $::bigtestdir/tmp/homwes $::bigtestdir/expected/yri_exome.homwes-out
} [subst {Files $::bigtestdir/tmp/homwes/homwes-out.work/homwes-out-NA19240.log and $::bigtestdir/expected/yri_exome.homwes-out/homwes-out.work/homwes-out-NA19240.log differ
child process exited abnormally}] error

cd $keepdir

testsummarize

