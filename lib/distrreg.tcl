proc distrreg_checkvalue {value} {
	if {$value in {{0} {}}} {
		return 0
	} elseif {$value in {chr chromosome 1 schr}} {
		return $value
	} elseif {[isint $value]} {
		return $value
	} elseif {[regexp {^s([0-9]+)$} $value]} {
		return $value
	} elseif {[file exists $value]} {
		return [file_absolute $value]
	} else {
		error "unknown value $value for -distrreg, must be an existing (region) file or one of: chr, chromosome, schr, schromosome, 1 or 0, a number or a number preceded by an s"
	}
}

proc distrreg_regs {regfile refseq} {
	if {$regfile eq "" || $regfile eq "0"} {
		return {}
	}
	set refseq [refseq $refseq]
	if {$regfile eq "1"} {
		set regfile chr
	}
	if {$regfile in "chr chromosome 1"} {
		set chromosomes [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai]
		unset -nocomplain a
		foreach chr $chromosomes {
			regsub {_.*$} $chr _ chr
			set a($chr) 1
		}
		return [bsort [array names a]]
	} elseif {$regfile in "schr schromosome"} {
		return [join [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai] " "]
	}
	if {[regexp {^s([0-9]+)$} $regfile temp regsize]} {
		set sequencedfile [gzfile [file dir $refseq]/extra/reg_*_sequencedgenome.tsv]
		if {![file exists $sequencedfile]} {
			error "sequencedfile not found ($sequencedfile)"
		}
		unset -nocomplain donea
		set result {}
		set list [split [string trim [cg select -sh /dev/null -f {chromosome begin end} $sequencedfile]] \n]
		foreach {cchr cbegin cend} [list_pop list 0] break
		set csize [expr {$cend-$cbegin}]
		list_foreach {chr begin end} $list {
			set size [expr {$end-$begin}]
			if {$chr ne $cchr || [expr {$csize + $size}] > $regsize} {
				if {[regexp {^[^_]*_} $cchr base]} {
					if {![info exists donea($base)]} {
						lappend result $base
						set donea($base) 1
					}
				} else {
					lappend result $cchr-$cbegin-$cend
				}
				set cchr $chr
				set cbegin $begin
				set cend $end
			} else {
				set cend $end
			}
			set csize [expr {$cend-$cbegin}]
		}
		lappend result $cchr-$cbegin-$cend
		return $result
	} elseif {[isint $regfile]} {
		set regsize $regfile
		unset -nocomplain donea
		set result {}
		list_foreach {chr size} [split [string trim [cg select -sh /dev/null -hp {chromosome size} -f {chromosome size} $refseq.fai]] \n] {
			if {[regexp {^[^_]*_} $chr base]} {
				if {![info exists donea($base)]} {
					lappend result $base
					set donea($base) 1
				}
				continue
			}
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
