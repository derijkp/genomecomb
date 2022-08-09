proc distrreg_checkvalue {value {default s50000000}} {
	if {$value in {{0} {}}} {
		return 0
	} elseif {$value eq "1"} {
		return $default
	} elseif {$value in {chr chromosome 1 schr}} {
		return $value
	} elseif {[isint $value]} {
		return $value
	} elseif {[regexp {^s([0-9]+)$} $value]} {
		return $value
	} elseif {[regexp {^r([0-9]+)$} $value]} {
		return $value
	} elseif {[regexp {^g([0-9]+)$} $value]} {
		return $value
	} elseif {[file exists $value]} {
		return [file_absolute $value]
	} else {
		error "unknown value $value for -distrreg, must be an existing (region) file or one of: chr, chromosome, schr, schromosome, 1 or 0, a number or a number preceded by an s, r, or g"
	}
}

proc distrreg_reg2loc {region {refseq {}}} {
	foreach {chr begin end} [split $region _-] break
	if {$begin eq ""} {set begin 0}
	if {$refseq eq ""} {
		if {$end eq ""} {set end 536870912}
	} elseif {$end eq ""} {
		set end [ref_chrsize $refseq $chr]
	}
	list $chr $begin $end
}

proc distrreg_reg2bed {bedfile regions {refseq {}}} {
	set o [open $bedfile w]
	foreach region $regions {
		foreach {chr begin end} [split $region :-] break
		if {$region eq "unaligned"} {
			continue
		}
		if {[isint $begin] && [isint $end]} {
			puts $o $chr\t$begin\t$end
			continue
		}
		global genomecomb_chrsizea
		if {![info exists genomecomb_chrsizea]} {
			list_foreach {tchr size} [split [string trim [cg select -sh /dev/null -hp {chromosome size} -f {chromosome size} $refseq.fai]] \n] {
				set genomecomb_chrsizea([chr_clip $tchr]) [list $tchr $size]
			}
		}
		set cchr [chr_clip $chr]
		if {[info exists genomecomb_chrsizea($cchr)]} {
			foreach {chr size} $genomecomb_chrsizea($cchr) break
			puts $o $chr\t0\t$size
		} elseif {[string index $chr end] eq "_"} {
			foreach tchr [array names genomecomb_chrsizea [chr_clip $chr]*] {
				foreach {chr size} $genomecomb_chrsizea($tchr) break
				puts $o $chr\t0\t$size
			}
		} else {
			puts $o $chr\t0\t536870912
		}
	}
	close $o
	return $bedfile
}

proc distrreg_reg2tsv {tsvfile regions {refseq {}}} {
	set o [open $tsvfile w]
	puts $o chromosome\tbegin\tend
	foreach region $regions {
		foreach {chr begin end} [split $region :-] break
		if {$region eq "unaligned"} {
			continue
		}
		if {[isint $begin] && [isint $end]} {
			puts $o $chr\t$begin\t$end
			continue
		}
		global genomecomb_chrsizea
		if {![info exists genomecomb_chrsizea]} {
			list_foreach {tchr size} [split [string trim [cg select -sh /dev/null -hp {chromosome size} -f {chromosome size} $refseq.fai]] \n] {
				set genomecomb_chrsizea([chr_clip $tchr]) [list $tchr $size]
			}
		}
		set cchr [chr_clip $chr]
		if {[info exists genomecomb_chrsizea($cchr)]} {
			foreach {chr size} $genomecomb_chrsizea($cchr) break
			puts $o $chr\t0\t$size
		} elseif {[string index $chr end] eq "_"} {
			foreach tchr [array names genomecomb_chrsizea [chr_clip $chr]*] {
				foreach {chr size} $genomecomb_chrsizea($tchr) break
				puts $o $chr\t0\t$size
			}
		} else {
			puts $o $chr\t0\t536870912
		}
	}
	close $o
	return $tsvfile
}

proc maxint args {
	set max [lindex $args 0]
	foreach num $args {
		if {$num > $max} {set max $num}
	}
	return $max
}

proc distrreg_use {distrreg {defdistrreg chr} {maxdistrreg {}}} {
	if {$distrreg == 0 || $maxdistrreg == 0} {
		return 0
	} elseif {$distrreg == 1} {
		return $defdistrreg
	} elseif {$distrreg eq "schr" || $maxdistrreg eq "schr"} {
		return schr
	} elseif {$distrreg eq "chr"} {
		return $distrreg
		if {$maxdistrreg eq "chr"} {
			return $maxdistrreg
		} else {
			return $distrreg
		}
	} elseif {[isint $distrreg]} {
		if {$maxdistrreg in "0 schr chr"} {
			return $maxdistrreg
		} elseif {[regexp {^s([0-9]+)$} $maxdistrreg temp maxnum]} {
			return s[maxint $maxnum $distrreg]
		} elseif {[isint $maxdistrreg]} {
			return [maxint $distrreg $maxdistrreg]
		} else {
			return $distrreg
		}
	} elseif {[regexp {^s([0-9]+)$} $distrreg temp num]} {
		if {$maxdistrreg in "0 schr chr"} {
			return $maxdistrreg
		} elseif {[regexp {^s([0-9]+)$} $maxdistrreg temp maxnum]} {
			return s[maxint $maxnum $num]
		} else {
			return $distrreg
		}
	} else {
		# filename
		if {$maxdistrreg ne ""} {
			return $maxdistrreg
		} else {
			return $distrreg
		}
	}
}

proc distrreg_addunaligned {regs {addunaligned 1}} {
	if {$addunaligned} {
		set last [lindex $regs end]
		set pre {}
		while {[bsort [list $last ${pre}unaligned]] ne [list $last ${pre}unaligned]} {
			append pre z
		}
		lappend regs ${pre}unaligned
	}
	return $regs
}

proc distrreg_regs {regfile refseq {addunaligned 1}} {
	if {$regfile eq "" || $regfile eq "0"} {
		return {}
	}
	set refseq [refseq $refseq]
	if {$regfile eq "1"} {
		set regfile chr
	}
	if {[regexp {^s([0-9]+)$} $regfile temp regsize]} {
		set sequencedfile [gzfile [file dir $refseq]/extra/reg_*_sequencedgenome.tsv]
		if {![file exists $sequencedfile]} {
			putslog "sequencedfile not found ($sequencedfile), using distrreg chr"
			set regfile chr
		} else {
			unset -nocomplain donea
			set result {}
			set list [split [string trim [cg select -sh /dev/null -f {chromosome begin end} $sequencedfile]] \n]
			foreach {cchr cbegin cend} [list_pop list 0] break
			set cbegin 0
			set cmax [ref_chrsize $refseq $cchr]
			set csize [expr {$cend-$cbegin}]
			list_foreach {chr begin end} $list {
				set size [expr {$end-$begin}]
				if {$chr ne $cchr} {
					if {[regexp {^[^_]*_} $cchr base]} {
						if {![info exists donea($base)]} {
							lappend result $base
							set donea($base) 1
						}
					} else {
						lappend result $cchr-$cbegin-$cmax
						set cmax [ref_chrsize $refseq $chr]
					}
					set cchr $chr
					set cbegin 0
					set cend $end
					set cmax [ref_chrsize $refseq $cchr]
				} elseif {[regexp {^[^_]*_} $cchr base]} {
					# do not add anything if chr1_*
				} elseif {[expr {$csize + $size}] > $regsize} {
					set mid [expr {($cend + $begin)/2}]
					lappend result $cchr-$cbegin-$mid
					set cchr $chr
					set cbegin $mid
					set cend $end
				} else {
					set cend $end
				}
				set csize [expr {$cend-$cbegin}]
			}
			lappend result $cchr-$cbegin-$cend
			set result [distrreg_addunaligned $result $addunaligned]
			return $result
		}
	}
	if {[regexp {^r([0-9]+)$} $regfile temp regsize]} {
		set ref [file tail [file dir $refseq]]
		set norep100000file [gzfile [file dir $refseq]/extra/reg_${ref}_norep100000.tsv]
		set fullgenomefile [gzfile [file dir $refseq]/extra/reg_*_fullgenome.tsv]
		if {![file exists $norep100000file]} {
			set rmskfile [gzfile [file dir $refseq]/reg_*_rmsk.tsv]
			if {![file exists $rmskfile]} {
				putslog "norep100000file $norep100000file not found, nor rmsk file ($rmskfile) to make it, using distrreg chr"
				set regfile chr
			} else {
				exec cg regjoin $rmskfile | cg select -f {{size=$end - $begin} *} \
					| cg select -q {$size > 100000} > $norep100000file.temp
				cg regsubtract $fullgenomefile $norep100000file.temp > $norep100000file
				file delete $norep100000file.temp
			}
		}
		unset -nocomplain donea
		set result {}
		set list [split [string trim [cg select -sh /dev/null -f {chromosome begin end} $norep100000file]] \n]
		foreach {cchr cbegin cend} [list_pop list 0] break
		set cbegin 0
		set cmax [ref_chrsize $refseq $cchr]
		set csize [expr {$cend-$cbegin}]
		list_foreach {chr begin end} $list {
			set size [expr {$end-$begin}]
			if {$chr ne $cchr} {
				if {[regexp {^[^_]*_} $cchr base]} {
					if {![info exists donea($base)]} {
						lappend result $base
						set donea($base) 1
					}
				} else {
					lappend result $cchr-$cbegin-$cmax
					set cmax [ref_chrsize $refseq $chr]
				}
				set cchr $chr
				set cbegin 0
				set cend $end
				set cmax [ref_chrsize $refseq $cchr]
			} elseif {[regexp {^[^_]*_} $cchr base]} {
				# do not add anything if chr1_*
			} elseif {[expr {$csize + $size}] > $regsize} {
				set mid [expr {($cend + $begin)/2}]
				lappend result $cchr-$cbegin-$mid
				set cchr $chr
				set cbegin $mid
				set cend $end
			} else {
				set cend $end
			}
			set csize [expr {$cend-$cbegin}]
		}
		lappend result $cchr-$cbegin-$cend
		set result [distrreg_addunaligned $result $addunaligned]
		return $result
	}
	if {[regexp {^g([0-9]+)$} $regfile temp regsize]} {
		set ref [file tail [file dir $refseq]]
		set nolowgenefile [gzfile [file dir $refseq]/extra/reg_${ref}_nolowgene250k.tsv]
		set fullgenomefile [gzfile [file dir $refseq]/extra/reg_*_fullgenome.tsv]
		if {![file exists $nolowgenefile]} {
			set temp [tempfile]
			cg regjoin $genedbtsv > $temp
			set temp2 [tempfile]
			cg regsubtract $fullgenomefile $temp > $temp2
			set temp3 [tempfile]
			cg select -overwrite 1 -q {$end - $begin > 200000} $temp2 $temp3
			# cg regsubtract $fullgenomefile $temp3 | cg select -f {{size=$end - $begin} *} | cg select -s size
			cg regsubtract $fullgenomefile $temp3 > $nolowgenefile.temp
			file rename $nolowgenefile.temp $nolowgenefile
		}
		unset -nocomplain donea
		set result {}
		set list [split [string trim [cg select -sh /dev/null -f {chromosome begin end} $nolowgenefile]] \n]
		foreach {cchr cbegin cend} [list_pop list 0] break
		set cbegin 0
		set cmax [ref_chrsize $refseq $cchr]
		set csize [expr {$cend-$cbegin}]
		list_foreach {chr begin end} $list {
			set size [expr {$end-$begin}]
			if {$chr ne $cchr} {
				if {[regexp {^[^_]*_} $cchr base]} {
					if {![info exists donea($base)]} {
						lappend result $base
						set donea($base) 1
					}
				} else {
					lappend result $cchr-$cbegin-$cmax
					set cmax [ref_chrsize $refseq $chr]
				}
				set cchr $chr
				set cbegin 0
				set cend $end
				set cmax [ref_chrsize $refseq $cchr]
			} elseif {[regexp {^[^_]*_} $cchr base]} {
				# do not add anything if chr1_*
			} elseif {[expr {$csize + $size}] > $regsize} {
				set mid [expr {($cend + $begin)/2}]
				lappend result $cchr-$cbegin-$mid
				set cchr $chr
				set cbegin $mid
				set cend $end
			} else {
				set cend $end
			}
			set csize [expr {$cend-$cbegin}]
		}
		lappend result $cchr-$cbegin-$cend
		set result [distrreg_addunaligned $result $addunaligned]
		return $result
	}
	if {$regfile in "chr chromosome 1"} {
		set chromosomes [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai]
		unset -nocomplain a
		foreach chr $chromosomes {
			regsub {_.*$} $chr _ chr
			set a($chr) 1
		}
		return [distrreg_addunaligned [bsort [array names a]] $addunaligned]
	} elseif {$regfile in "schr schromosome"} {
		return [distrreg_addunaligned [join [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai] " "] $addunaligned]
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
		set result [distrreg_addunaligned $result $addunaligned]
		return $result
	}
	set f [gzopen $regfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set result {}
	while {[gets $f line] != -1} {
		lappend result [join [list_sub [split $line \t] $poss] -]
	}
	set result [distrreg_addunaligned $result $addunaligned]
	return $result
}
