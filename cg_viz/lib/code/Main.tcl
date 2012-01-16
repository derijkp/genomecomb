proc main args {
#putsvars ::appdir
lappend ::auto_path [file normalize $::Classy::appdir/../lib]
mainw .mainw
focus .mainw
#Classy::cmd
if {[file exists [lindex $args 0]]} {
	.mainw opentsv {*}$args
} else {
	# .mainw opendb {*}$args
}
}


if 0 {
	cd ~/dev/completegenomics/tests
	lappend auto_path ~/dev/completegenomics/cg_viz/lib  ~/dev/completegenomics/cg_viz/lib-exp ~/dev/completegenomics/cg_viz/lib/code ~/dev/completegenomics/cg_viz/lib/interface
	file copy data/expected-annotate-vars_annottest-gene_test.tsv test.tsv
	set file test.tsv
	set object .mainw.tb
	set argv [list $file]
	source ../cg_viz/cg_viz.tcl
	$object open $file

set file /complgen/projects/try/temp2.annot.temp

	set object tb
	table_tsv new tb
	tb open $file
}
