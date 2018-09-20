proc distrreg_regs {regfile refseq} {
	if {$regfile in "chr chromosome 1"} {
		return [join [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai] " "]
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
