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
	} elseif {$value eq "g"} {
		return $value
	} elseif {[regexp {^g([0-9]+)$} $value]} {
		return $value
	} elseif {[regexp {^x([0-9]+)$} $value]} {
		return $value
	} elseif {[file exists $value]} {
		return [file_absolute $value]
	} else {
		error "unknown value $value for -distrreg, must be an existing (region) file or one of: chr, chromosome, schr, schromosome, 1 or 0, g, a number or a number preceded by an s, r, or g"
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
	} elseif {[regexp {^x([0-9]+)$} $distrreg temp num]} {
		if {$maxdistrreg in "0 schr chr"} {
			return $maxdistrreg
		} elseif {[regexp {^([sgrx])([0-9]+)$} $maxdistrreg temp maxtype maxnum]} {
			return $maxtype[maxint $maxnum $num]
		} else {
			return $distrreg
		}
	} elseif {[regexp {^([sgrx])([0-9]+)$} $distrreg temp type num]} {
		if {$maxdistrreg in "0 schr chr"} {
			return $maxdistrreg
		} elseif {[regexp {^([sgrx])([0-9]+)$} $maxdistrreg temp maxtype maxnum]} {
			return $type[maxint $maxnum $num]
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

proc distrreg_norep100000file {refdir} {
	set ref [file tail $refdir]
	set norep100000file [gzfile $refdir/extra/reg_${ref}_norep100000.tsv.zst]
	if {![file exists $norep100000file]} {
		set fullgenomefile [gzfile $refdir/extra/reg_*_fullgenome.tsv]
		set rmskfile [gzfile $refdir/reg_*_rmsk.tsv]
		if {![file exists $rmskfile]} {
			putslog "norep100000file $norep100000file not found, nor rmsk file ($rmskfile) to make it, using distrreg chr"
		} else {
			exec cg regjoin $rmskfile | cg select -f {{size=$end - $begin} *} \
				| cg select -q {$size > 100000} > $norep100000file.temp
			cg regsubtract $fullgenomefile $norep100000file.temp | cg zst > $norep100000file.temp2
			file delete $norep100000file.temp
			file rename $norep100000file.temp2 $norep100000file
			cg zstindex $norep100000file
		}
	}
	return $norep100000file
}

proc distrreg_nolowgene {refdir {cutoff 200000} {name {}}} {
	if {$cutoff eq "" && $name eq "" && [file exists $refdir/extra/reg_${ref}_nolowgene.tsv.zst]} {
		return $refdir/extra/reg_${ref}_nolowgene.tsv.zst
	}
	if {$cutoff eq ""} {
		set cutoff 200000
	}
	if {$name eq ""} {
		set name $cutoff
		regsub {000$} $name k name
	}
	set ref [file tail $refdir]
	set nolowgenefile [gzfile $refdir/extra/reg_${ref}_nolowgene$name.tsv.zst]
	if {![file exists $nolowgenefile]} {
		set temp [lindex [bsort [gzfiles $refdir/extra/reg_${ref}_nolowgene*.tsv.zst]] 0] 
		if {[file exists $temp]} {
			return $temp
		}
	}
	if {![file exists $nolowgenefile]} {
		set fullgenomefile [gzfile $refdir/extra/reg_*_fullgenome.tsv]
		set temp [tempfile]
		set genedbtsv [ref_tsvtranscripts $refdir]
		cg regjoin $genedbtsv > $temp
		set temp2 [tempfile]
		cg regsubtract $fullgenomefile $temp > $temp2
		# cg select -f {chromosome {size=$end - $begin}} $temp2 | cg select -s -size | head -40
		set temp3 [tempfile]
		cg select -overwrite 1 -q "\$end - \$begin > $cutoff" $temp2 $temp3
		# cg regsubtract $fullgenomefile $temp3 | cg select -f {{size=$end - $begin} *} | cg select -s size
		cg regsubtract $fullgenomefile $temp3 | cg zst > $nolowgenefile.temp
		file rename $nolowgenefile.temp $nolowgenefile
		cg zstindex $nolowgenefile
	}
	return $nolowgenefile
}

proc cg_distrreg_nolowgene {args} {
	distrreg_nolowgene {*}$args
}

proc distrreg_group_read {refseq groupchraVar {elementsaVar {}}} {
	upvar $groupchraVar groupchra
	if {$elementsaVar ne ""} {upvar $elementsaVar elementsa}
	unset -nocomplain groupchra
	unset -nocomplain elementsa
	if {[file exists $refseq.groupchromosomes]} {
		set f [open $refseq.groupchromosomes]
		set header [tsv_open $f]
		while {[gets $f line] != -1} {
			foreach {chr group} [split $line \t] break
			set groupchra($chr) $group
			lappend elementsa($group) $chr
		}
	}
}

proc distrreg_group_base {groupchraVar chr} {
	upvar $groupchraVar groupchra
	
	if {[info exists groupchra($chr)]} {
		if {$groupchra($chr) ne $chr} {
			set base $groupchra($chr)
		} else {
			set base {}
		}
	} elseif {[regexp {^[^_]*_} $chr base]} {
	} else {
		set base {}
	}
	return $base
}

proc distrreg_regs {regfile refseq {type s} {addunaligned 1}} {
	if {$regfile eq "" || $regfile eq "0"} {
		return {}
	}
	set refseq [refseq $refseq]
	if {$regfile eq "1"} {
		set regfile chr
	}
	if {[regexp {^x([0-9]+)$} $regfile temp regsize]} {
		set regfile $type$regsize
	}
	distrreg_group_read [refseq $refseq] groupchra
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
				set base [distrreg_group_base groupchra $cchr]
				if {$chr ne $cchr} {
					if {$base ne ""} {
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
				} elseif {$base ne ""} {
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
		set norep100000file [distrreg_norep100000file [file dir $refseq]]
		if {![file exists $norep100000file]} {
			putslog "norep100000file not found ($norep100000file), using distrreg chr"
			set regfile chr
		} else {
			unset -nocomplain donea
			set result {}
			set list [split [string trim [cg select -sh /dev/null -f {chromosome begin end} $norep100000file]] \n]
			foreach {cchr cbegin cend} [list_pop list 0] break
			set cbegin 0
			set cmax [ref_chrsize $refseq $cchr]
			set csize [expr {$cend-$cbegin}]
			list_foreach {chr begin end} $list {
				set size [expr {$end-$begin}]
				set base [distrreg_group_base groupchra $cchr]
				if {$chr ne $cchr} {
					if {$base ne ""} {
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
				} elseif {$base ne ""} {
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
	if {$regfile eq "g"} {
		set distrgfile [gzfile [file dir $refseq]/extra/reg_*_distrg.tsv]
		if {[file exists $distrgfile]} {
			set f [open $distrgfile]
			set header [tsv_open $f]
			set result {}
			while {[gets $f line] != -1} {
				foreach {chr b e} {{} {} {}} break
				foreach {chr b e} [split $line \t] break
				if {$b eq ""} {
					lappend result $chr
				} else {
					lappend result $chr-$b-$e
				}
			}
			set result [distrreg_addunaligned $result $addunaligned]
			return $result
		}
		set regfile g5000000
	}
	if {[regexp {^g([0-9]+)$} $regfile temp regsize]} {
		set nolowgenefile [distrreg_nolowgene [file dir $refseq]]
		if {![file exists $nolowgenefile]} {
			putslog "nolowgenefile not found ($nolowgenefile), using distrreg chr"
			set regfile chr
		} else {
			unset -nocomplain donea
			set result {}
			set list [split [string trim [cg select -sh /dev/null -f {chromosome begin end} $nolowgenefile]] \n]
			foreach {cchr cbegin cend} [list_pop list 0] break
			set cbegin 0
			set cmax [ref_chrsize $refseq $cchr]
			set csize [expr {$cend-$cbegin}]
			list_foreach {chr begin end} $list {
				set size [expr {$end-$begin}]
				set base [distrreg_group_base groupchra $cchr]
				if {$chr ne $cchr} {
					if {$base ne ""} {
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
				} elseif {$base ne ""} {
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
	if {$regfile in "chr chromosome 1"} {
		set chromosomes [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai]
		unset -nocomplain a
		foreach chr $chromosomes {
			set base [distrreg_group_base groupchra $chr]
			if {$base ne ""} {
				set a($base) 1
			} else {
				set a($chr) 1
			}
		}
		return [distrreg_addunaligned [bsort [array names a]] $addunaligned]
	} elseif {$regfile in "schr schromosome"} {
		return [distrreg_addunaligned [join [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai] " "] $addunaligned]
	} elseif {[isint $regfile]} {
		set regsize $regfile
		unset -nocomplain donea
		set result {}
		list_foreach {chr size} [split [string trim [cg select -sh /dev/null -hp {chromosome size} -f {chromosome size} $refseq.fai]] \n] {
			set base [distrreg_group_base groupchra $chr]
			if {$base ne ""} {
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
