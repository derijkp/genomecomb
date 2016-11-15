#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set basedir /data/genomecomb.testdata

set keepdir [pwd]

test homwes {} {
	test_cleantmp
	set varfile $::basedir/exomes_yri/compar/annot_compar-yri_exome_test.tsv.lz4
	set dbdir $::basedir/refseq/hg19_test
	cg homwes -v 1 -dbdir $dbdir $varfile NA19240 tmp/homwes-out.tsv 2> /dev/null
	catch {exec diff -rq tmp $::basedir/expected/yri_exome.homwes-out} msg
	set msg
} [subst {Files tmp/homwes-out.work/homwes-out-NA19240.log and $::basedir/expected/yri_exome.homwes-out/homwes-out.work/homwes-out-NA19240.log differ
child process exited abnormally}]

cd $keepdir

testsummarize

