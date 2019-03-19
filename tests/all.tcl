#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" "$@"

proc runtests file {
	cd $::appdir/tests
	uplevel source $file
}

runtests analysis.tcl
runtests annot.tcl
runtests bcol.tcl
runtests cg2tsv.tcl
runtests clip.tcl
runtests convert.tcl
runtests compress.tcl
runtests fastq.tcl
runtests genome_seq.tcl
runtests primercheck.tcl
runtests libext.tcl
runtests lift.tcl
runtests mirannot.tcl
runtests misc.tcl
runtests reg.tcl
runtests seq.tcl
runtests remap.tcl
runtests select.tcl
runtests select_group.tcl
runtests select_saggregates.tcl
runtests multiselect.tcl
runtests tsv.tcl
runtests val.tcl
runtests varia.tcl
runtests vcf.tcl
runtests vcf2tsv.tcl
runtests tsv2vcf.tcl
runtests map.tcl
runtests var.tcl
runtests job.tcl
runtests compar.tcl
runtests nanopore.tcl
runtests svmulticompar.tcl
runtests sv.tcl
runtests pmulticompar.tcl
runtests process_multicompar.tcl

# next ones take longer, use larger data in genomecomb.testdata
runtests queries.tcl
runtests reports.tcl
runtests bam.tcl
runtests homwes.tcl

# long time, run separately
# runtests process.tcl

# not really used
# runtests mselect.tcl
