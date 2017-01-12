#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]

test homwes {} {
	file delete -force $::bigtestdir/tmp/homwes
	file mkdir $::bigtestdir/tmp/homwes
	set varfile $::bigtestdir/ori/exomes_yri/compar/annot_compar-yri_exome_test.tsv.lz4
	set dbdir $::bigtestdir/refseq/hg19_test
	cg homwes --stack 1 -v 1 -dbdir $dbdir $varfile NA19240 $::bigtestdir/tmp/homwes/homwes-out.tsv 2> /dev/null
	catch {exec diff -rq $::bigtestdir/tmp/homwes $::bigtestdir/expected/yri_exome.homwes-out} msg
	set msg
} [subst {Files $::bigtestdir/tmp/homwes/homwes-out.work/homwes-out-NA19240.log and $::bigtestdir/expected/yri_exome.homwes-out/homwes-out.work/homwes-out-NA19240.log differ
child process exited abnormally}]

cd $keepdir

testsummarize

