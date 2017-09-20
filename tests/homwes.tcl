#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

set keepdir [pwd]

test homwes {basic} {
	file delete -force $::bigtestdir/tmp/homwes
	file mkdir $::bigtestdir/tmp/homwes
	set varfile $::bigtestdir/ori/exomes_yri.ori/compar/annot_compar-yri_exome_test.tsv.lz4
	set dbdir $::bigtestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile NA19240 $::bigtestdir/tmp/homwes/homwes-out.tsv
	cg tsvdiff -q 1 -x *.log $::bigtestdir/tmp/homwes $::bigtestdir/expected/yri_exome.homwes-out
} {}

test homwes {sample not given -> multiple samples} {
	file delete -force $::bigtestdir/tmp/homwes
	file mkdir $::bigtestdir/tmp/homwes
	set varfile $::bigtestdir/ori/exomes_yri.ori/compar/annot_compar-yri_exome_test.tsv.lz4
	set dbdir $::bigtestdir/refseqtest/hg19
	cg homwes --stack 1 -dbdir $dbdir $varfile {} $::bigtestdir/tmp/homwes_multi/homwes-out.tsv
	cg tsvdiff -q 1 -x *.log -x log_jobs $::bigtestdir/tmp/homwes_multi $::bigtestdir/expected/homwes_multi
} {}

cd $keepdir

testsummarize

