package require Extral
catch {tk appname test}

set bigtestdir /data/genomecomb.testdata

package require pkgtools
namespace import -force pkgtools::*
package require Extral

proc test {args} {
	test_cleantmp
	pkgtools::test {*}$args
	cd $::testdir
	return {}
}

proc checkdiff args {
	global e
	set err [catch {exec diff {*}$args} e]
	if {$err && $e ne {child process exited abnormally}} {error $e}
}

# pkgtools::testleak 100

set keeppath $::env(PATH)
set script [info script] ; if {$script eq ""} {set script ./t}
set testdir [file dir [file normalize $script]]
set refseqdir $testdir/genomecomb.testdata/refseq
set appdir [file dir [file dir [file normalize $script]]]
append ::env(PATH) :$appdir/bin
# putsvars ::env(PATH)
set env(SCRATCHDIR) [file dir [tempdir]]
source $appdir/lib/file.tcl ; pwd

proc test_cleantmp {} {
	foreach file [glob -nocomplain $::testdir/tmp/*] {
		catch {file attributes $file -permissions ugo+xw}
		catch {file delete -force $file}
	}
	cg indexclean
}

proc write_tab {file data {comment {}}} {
	set data [split [string trim $data] \n]
	set f [open $file w]
	if {$comment ne ""} {puts $f $comment}
	foreach line $data {
		puts $f [join $line \t]
	}
	close $f
}

proc diff_tab {file data {comment {}}} {
	set tempfile [tempfile]
	set data [split [string trim $data] \n]
	set f [open $tempfile w]
	if {$comment ne ""} {puts $f $comment}
	foreach line $data {
		puts $f [join $line \t]
	}
	close $f
}

proc test_genomecombdir {} {
	set expdir tmp/test
	file mkdir $expdir/samples
	file mkdir $expdir/compar
	set samplesdir $expdir/samples
	foreach n {1 2 3} {
		file mkdir $samplesdir/sample$n
	}
	cg splitalleles data/var_annot.sft > $samplesdir/sample1/var-sample1.tsv
	cg splitalleles data/var_annot2.sft > $samplesdir/sample2/var-sample2.tsv
	cg splitalleles data/var_annot2seq.sft > $samplesdir/sample3/prevar-sample3.tsv
	cg select -f {sequenced *} $samplesdir/sample3/prevar-sample3.tsv $samplesdir/sample3/var-sample3.tsv
	file copy data/sreg-annot1.sft $samplesdir/sample1/sreg-sample1.tsv
	file copy data/sreg-annot2.sft $samplesdir/sample2/sreg-sample2.tsv
	file copy data/sreg-annot2.sft $samplesdir/sample3/sreg-sample3.tsv
	cg multicompar -reannot -split 1 $expdir/compar/compar-test.tsv {*}[lsort -dict [glob $samplesdir/*/var-*.tsv]]
	file delete $expdir/compar/compar-test.tsv.reannot
	file delete $expdir/compar/compar-test.tsv.old
	# exec diff $expdir/compar/compar-test.tsv data/expected-multicompar-split-reannot.sft
}

lappend auto_path $appdir/lib $appdir/lib-exp $appdir/libext

# remove tmp if it is a unexisting link
if {![file exists tmp]} {catch {file delete tmp}}
file mkdir tmp
test_cleantmp
set dbopt {}
set dboptt {}
