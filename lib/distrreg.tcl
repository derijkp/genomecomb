proc distrreg_regs {regfile refseq} {
	if {$regfile in "chr chromosome 1"} {
		return [exec cut -f 1 $refseq.fai]
	}
	set f [gzopen $regfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set result {}
	while {[gets $f line] != -1} {
		lappend result [join [list_sub [split $line \t] $poss] _]
	}
	return $result
}
