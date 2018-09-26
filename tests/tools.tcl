package require Extral
catch {tk appname test}

set bigtestdir /data/genomecomb.testdata
set smalltestdir /data/genomecomb.smalltestdata
if {![info exists refseqdir]} {
	set refseqdir /data/genomecomb.smalltestdata/refseqtest
}

package require pkgtools
namespace import -force pkgtools::*
package require Extral

set test_cleantmp 1

proc test {args} {
	if {[get ::test_cleantmp 1]} {test_cleantmp}
	catch {job_init}
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
set appdir [file dir [file dir [file normalize $script]]]
source $appdir/lib/file.tcl
append ::env(PATH) [pathsep]$appdir/bin
# putsvars ::env(PATH)
set env(SCRATCHDIR) [file dir [tempdir]]
source $appdir/lib/file.tcl ; pwd

proc test_cleantmp {} {
	foreach file [list_remove [glob -nocomplain $::testdir/tmp/* $::testdir/tmp/.*] $::testdir/tmp/.. $::testdir/tmp/.] {
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

proc file_regsub {exp subSpec file resultfile} {
	set f [open $file]
	set o [open $resultfile w]
	while {[gets $f line] != -1} {
		regsub -all $exp $line $subSpec line
		puts $o $line
	}
	close $o
	close $f
}

proc write_sam {file data} {
	# creates a sam file based on data
	# sequences are all one base (given in data or A) and qualities are all -
	set tempfile [tempfile]
	set o [open $tempfile w]
	set num 1
	foreach line [split [string trim $data] \n] {
		set base A
		foreach {chr1 pos1 cigar1 seq1 chr2 pos2 cigar2 seq2 base} $line break
		if {$base eq ""} {set base A}
		if {[isint $seq1]} {
			set size1 $seq1
			set seq1 [string_fill $base $size1]
		} else {
			set size1 [string length $seq1]
		}
		set qual1 [string_fill - $size1]
		if {[isint $seq2]} {
			set size2 $seq2
			set seq2 [string_fill $base $size2]
		} else {
			set size2 [string length $seq2]
		}
		set qual2 [string_fill - $size2]
		if {$chr2 eq $chr1} {set c2 =} else {set c2 $chr2}
		puts $o [join [list A$num 99 $chr1 $pos1 60 $cigar1 $c2 $pos2 [expr {$pos2+$size2-$pos1}] $seq1 $qual1 RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25] \t]
		if {$chr2 eq $chr1} {set c1 =} else {set c1 $chr1}
		puts $o [join [list A$num 147 $chr2 $pos2 60 $cigar2 $c1 $pos1 -[expr {$pos2+$size2-$pos1}] $seq2 $qual2 RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25] \t]
		incr num
	}
	close $o
	set o [open $file w]
	foreach line [split [string trim {
		@HD	VN:1.4	GO:none	SO:coordinate
		@SQ	SN:chr1	LN:249250621
		@SQ	SN:chr2	LN:243199373
		@SQ	SN:chr3	LN:198022430
		@SQ	SN:chr4	LN:191154276
		@SQ	SN:chr5	LN:180915260
		@SQ	SN:chr6	LN:171115067
		@SQ	SN:chr7	LN:159138663
		@SQ	SN:chr8	LN:146364022
		@SQ	SN:chr9	LN:141213431
		@SQ	SN:chr10	LN:135534747
		@SQ	SN:chr11	LN:135006516
		@SQ	SN:chr12	LN:133851895
		@SQ	SN:chr13	LN:115169878
		@SQ	SN:chr14	LN:107349540
		@SQ	SN:chr15	LN:102531392
		@SQ	SN:chr16	LN:90354753
		@SQ	SN:chr17	LN:81195210
		@SQ	SN:chr18	LN:78077248
		@SQ	SN:chr19	LN:59128983
		@SQ	SN:chr20	LN:63025520
		@SQ	SN:chr21	LN:48129895
		@SQ	SN:chr22	LN:51304566
		@SQ	SN:chrM	LN:16571
		@SQ	SN:chrX	LN:155270560
		@SQ	SN:chrY	LN:59373566
		@RG	ID:sample1	PL:illumina	PU:sample1	LB:solexa-123	SM:sample1
		@PG	ID:GATK IndelRealigner	VN:2.4-9-g532efad	{CL:knownAlleles=[] targetIntervals=test.intervals LODThresholdForCleaning=5.0 consensusDeterminationModel=USE_READS entropyThreshold=0.15 maxReadsInMemory=150000 maxIsizeForMovement=3000 maxPositionalMoveAllowed=200 maxConsensuses=30 maxReadsForConsensuses=120 maxReadsForRealignment=20000 noOriginalAlignmentTags=false nWayOut=null generate_nWayOut_md5s=false check_early=false noPGTag=false keepPGTags=false indelsFileForDebugging=null statisticsFileForDebugging=null SNPsFileForDebugging=null}
	}] \n] {
		puts $o [join $line \t]
	}
	close $o
	exec gnusort8 -t \t -N -s -k3,3 -k4,4 -k1,1 $tempfile >> $file
}

proc diffanalysisinfo {file1 file2} {
	set file1 [lindex [glob $file1] 0]
	set file2 [lindex [glob $file2] 0]
	catch {exec cg tsvdiff -t xl $file1 $file2 | grep -v _version} temp
	set len [llength [split $temp \n]]
	if {$len != 3 && $len != 1} {error "error comparing $file1 and $file2: $temp"}
}

lappend auto_path $appdir/lib $appdir/lib-exp $appdir/libext

# remove tmp if it is a unexisting link
if {![file exists tmp]} {catch {file delete tmp}}
file mkdir tmp
set dbopt {}
set dboptt {}
