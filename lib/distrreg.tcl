proc distrreg_regs {regfile refseq} {
	if {$regfile in "chr chromosome 1"} {
		return [join [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai] " "]
	}
	if {[isint $regfile]} {
		set regsize $regfile
		set result {}
		list_foreach {chr size} [split [string trim [cg select -sh /dev/null -hp {chromosome size} -f {chromosome size} $refseq.fai]] \n] {
			for {set pos 0} {$pos < $size} {incr pos $regsize} {
				set end [expr {$pos+$regsize}]
				if {$end > $size} {set end $size}
				lappend result $chr-$pos-$end
			}
		}
		return $result
	}
	set f [gzopen $regfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set result {}
	while {[gets $f line] != -1} {
		lappend result [join [list_sub [split $line \t] $poss] -]
	}
	return $result
}
