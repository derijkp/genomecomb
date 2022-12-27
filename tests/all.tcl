#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" "$@"

set script [info script] ; if {$script eq ""} {set script ./t}
set appdir [file dir [file dir [file normalize $script]]]
# set script [pwd]/all.tcl ; set appdir [file dir [pwd]] ; set argv {}
puts "running all tests with script=$script appdir=$appdir argv=$argv"
if {[inlist $argv testsge]} {
	puts stderr "testing sge"
	set testsge 1
} else {
	set testsge 0
}
set keeptimes {}
proc runtests file {
	cd $::appdir/tests
	set time [time {uplevel source $file}]
	set time [expr {[lindex $time 0]/1000000.0}]
	lappend ::keeptimes $file\t$time
	puts "ran $file: $time seconds"
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
runtests select_aaggregates.tcl
runtests multiselect.tcl
runtests tsv.tcl
runtests index.tcl
runtests varia.tcl
runtests bsort.tcl
runtests val.tcl
runtests vcf.tcl
runtests vcf2tsv.tcl
runtests tsv2vcf.tcl
runtests compar.tcl
runtests rna.tcl
runtests svmulticompar.tcl
runtests nanopore.tcl
runtests queries.tcl
runtests reports.tcl
runtests bam.tcl
runtests homwes.tcl
runtests pmulticompar.tcl
runtests process_multicompar.tcl
runtests process_local.tcl
runtests bam_clean.tcl
runtests job.tcl
runtests map.tcl
runtests sv.tcl
runtests var.tcl
# not yet updated
# runtests meth.tcl

puts "all tests finished"
puts "times (seconds):\n[join $::keeptimes \n]"
testsummarize

# take long time, run separately
# ./process_small.tcl
# ./process_ont.tcl

# take longer still, run separately on cluster (will be run in ~/genomecomb_giab_testdata)
# without parameter, the code will only check previous runs (should be run after analysis is finished on the cluster)
# use with parameter run to delete results in tmp and completely (re)run the analysis
#./process_giab.tcl run
#./process_giab.tcl
#./process_giab_ont.tcl run
#./process_giab_ont.tcl
# not yet working adapted to this way of working:
#./process_giab_exome.tcl run
#./process_giab_exome.tcl

# larger still, not kept up to date
# ./process.tcl
# ./process_large.tcl
# not really used
# runtests mselect.tcl
