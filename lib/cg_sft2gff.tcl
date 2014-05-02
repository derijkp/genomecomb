proc cg_sft2gff {args} {
	if {[llength $args] > 2} {
		puts stderr "format is cg sft2gff ?tsvfile? ?gfffile?"
		# errorformat select
		exit 1
	}
	foreach {tsvfile gfffile} $args break
	if {$tsvfile eq ""} {
		set f stdin
	} else {
		set f [gzopen $tsvfile]
	}
	if {$gfffile eq ""} {
		set o stdout
	} else {
		set o [open $gfffile w]
	}
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 3]
	set poss [linsert $poss 1 [lindex [list_cor $header {source}] 0]]
	set poss [linsert $poss 2 [lindex [list_cor $header {type}] 0]]
	lappend poss [lindex [list_cor $header {score}] 0]
	lappend poss [lindex [list_cor $header {strand dir}] 0]
	lappend poss [lindex [list_cor $header {phase}] 0]
	# puts -nonewline $o $comment
	# puts $o "# converted to gff from $tsvfile"
	set transmap [list \; %3B = %3D % %25 & %26 , %2C \t %09 \n %0A]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set result [list_sub $line $poss]
		if {[lindex $result 6] eq ""} {lset result 6 +}
		lset result 0 [string map $transmap [lindex $result 0]]
		lset result 3 [expr {[lindex $result 3] + 1}]
		lset result 4 [expr {[lindex $result 4] + 1}]
		lappend result {}
		set result [list_change $result {{} .}]
		puts $o [join $result \t]
	}
	if {$o ne "stdout"} {close $o}
	if {$f ne "stdout"} {close $f}
}

if 0 {
	set tsvfile tests/data/reg4.tsv
	set tsvfile tests/data/testvars.tsv
	set o stdout
}
