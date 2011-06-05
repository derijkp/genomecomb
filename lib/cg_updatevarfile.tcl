#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

if 0 {
	cd /complgen/projects/dlb1
	set file oldcmt71_compar.tsv.old
}

proc cg_updatevarfile {args} {
	if {([llength $args] != 3)} {
		errorformat updatevarfile
		exit 1
	}
	foreach {file resultfile dbdir} $args break
	if {[file exists $resultfile]} {error "$resultfile exists"}
	catch {close $o} ; catch	{close $f}
	set f [open $file]
	set o [open $resultfile.temp w]
	set header [split [gets $f] \t]
	set nheader [list_concat {chromosome begin end type ref alt} [list_remove $header chromosome begin end type ref alt]]
	set fposs [tsv_basicfields $header 3]
	set poss [list_cor $header $nheader]
	set poss [lreplace $poss 0 2 {*}$fposs]
	set aposs [list_find -glob $nheader alleleSeq*]
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	puts $o [join $nheader \t]
	set count 0
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set line [list_sub $line $poss]
		foreach {chr start end type ref alt} $line break
		incr count
		if {$count > 1000000} {
			putslog $chr:$start-$end
			set count 0
		}
		set gref [string toupper [genome_get $fg $chr $start $end]]
		if {$gref ne $ref && $ref ne ""} {error "different ref for line:\n$line"}
		lset line 4 $gref
		if {$doalt} {
			set alleles [list_sub $line $aposs]
			set alt [list_remove [list_remdup $alleles] - ? N $gref {}]
			if {[llength $alt] == 0} {set alt ?}
			lset line 5 [join $alt ,]
		}
		puts $o [join $line \t]
	}
	close $o
	close $f
	file rename -force $resultfile.temp $resultfile
}


