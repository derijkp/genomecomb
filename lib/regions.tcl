#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc open_region {f {headerVar {}}} {
	if {[string length $headerVar]} {
		upvar $headerVar header
	}
	set header [tsv_open $f]
	if {[string index $header 0] eq "#"} {
		set header [string range $header 1 end]
	}
	set poss2 [tsv_basicfields $header 3]
	return $poss2
}

proc get_region {f poss} {
	while 1 {
		set line [split [gets $f] \t]
		if {[llength $line]} break
		if {[eof $f]} break
	}
	set result [list_sub $line $poss]
}

proc refconsregions {varfile} {
	putslog "getting ref-(in)consistent regions from $varfile"
	cg select -q {$varType == "ref-consistent" || $varType == "ref-inconsistent" || $varType == "no-call-rc" || $varType == "no-call-ri"} $varfile rctemp.tsv
	cg regjoin rctemp.tsv >@stdout
	file delete rctemp.tsv
}

proc cg_refconsregions {args} {
	global scriptname action
	if {[llength $args] != 1} {
		error "format is: $scriptname $action variation_file\n - outputs regions annotated as ref-(in)consistent from variation file"
	}
	foreach {varfile} $args break
	refconsregions $varfile
}

proc nocallregions {varfile outfile} {
	putslog "getting partial no-call regions from $varfile"
	set h [cg select -h $varfile]
	if {[inlist $h allele]} {
		cg select -q {$varType == "no-call" && $allele != "all"} $varfile nctemp.tsv
	} else {
		cg select -q {$varType == "no-call" && $haplotype != "all"} $varfile nctemp.tsv
	}
	cg regjoin nctemp.tsv > $outfile
	file delete nctemp.tsv
}

proc findregionfile {file} {
	set tail [file tail $file]
	regsub {^[^-]*-} $tail {} tail
	return [file dir $file]/sreg-$tail
}

proc regions2bed {regions refseq} {
	set result {}
	foreach region $regions {
		foreach r [samregions $region $refseq 1] {
			foreach {c b e} [split $r :-] break
			incr b -1
			lappend $c\t%b\t$e
		}
	}
	return [join $result \n]\n
}

proc samregions {region {refseq {}} {full 0}} {
	if {$region eq ""} {return $region}
	set split [split $region :-]
	foreach {chr begin end} {{} {} {}} break
	foreach {chr begin end} $split break
	if {$refseq ne ""} {
		distrreg_group_read [refseq $refseq] groupchra elementsa
	}
	if {!$full && ($begin eq "" || $end eq "")} {
		if {$begin ne "" || $end ne ""} {
			error "incorrect region:, must be either chr or chr:begin-end"
		}
		if {[info exists elementsa($chr)]} {
			return $elementsa($chr)
		} elseif {[regexp _$ $chr]} {
			set refseq [refseq $refseq]
			set chromosomes [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai]
			set result {}
			foreach tchr $chromosomes {
				if {[regexp ^$chr $tchr]} {
					lappend result $tchr
				}
			}
			return $result
		} else {
			return [list $chr]
		}
	}
	if {[info exists elementsa($chr)]} {
		set chromosomes $elementsa($chr)
		set result {}
		foreach tchr $chromosomes {
			lappend result $tchr:1-[ref_chrsize $refseq $tchr]
		}
		return $result
	} elseif {[regexp _$ $chr]} {
		set chromosomes [cg select -sh /dev/null -hp {chromosome size p1 p2} -f chromosome $refseq.fai]
		set result {}
		foreach tchr $chromosomes {
			if {[regexp ^$chr $tchr]} {
				lappend result $tchr:1-[ref_chrsize $refseq $tchr]
			}
		}
		return $result
	}
	#if {$begin eq ""} {set begin 1} else {incr begin}
	if {$begin eq ""} {set begin 1}
	if {$end eq ""} {
		set refseq [refseq $refseq]
		set end [ref_chrsize $refseq $region]
	}
	return [list $chr:$begin-$end]
}

proc samregions_gatk {region {refseq {}} {full 0}} {
	set regions [samregions $region $refseq $full]
	if {[llength $regions] == 1} {return [lindex $regions 0]}
	set tempbed [tempfile].bed
	distrreg_reg2bed $tempbed $regions $refseq
	return $tempbed
}

proc samregion {region {refseq {}} {full 0}} {
	set regions [samregions $region $refseq $full]
	if {[llength $regions] != 1} {
		error "error getting samregion from $region: matches [llength $regions] regions iso 1"
	}
	lindex $regions 0
}

proc regions_insert_next {list posVar cVar bVar eVar} {
	upvar $posVar pos
	upvar $cVar c
	upvar $bVar b
	upvar $eVar e
	set cur [lindex $list $pos]
	foreach {c b e} {{} {} {}} break
	foreach {c b e} [split $cur :-] break
	incr pos
	return $cur	
}

proc regions_insert {regions rDNA refseq} {
	set regions [bsort $regions]
	set inserts {}
	foreach r [bsort $rDNA] {
		lappend inserts {*}[samregions $r $refseq 1]
	}
	set inserts [bsort $inserts]
	set pos 0
	set cur [regions_insert_next $inserts pos curc curb cure]
	set regionpos 0
	set region [regions_insert_next $regions regionpos c b e]
	set result {}
	while 1 {
		if {$cur eq "" && $region eq ""} break
		if {$c eq $curc} {
			if {$b eq ""} {
				# full chromosome overlap
				foreach {c b e} [split [samregions $region $refseq 1] :-] break
				if {$b < $curb} {
					lappend result $c:$b-$curb
				}
				lappend result $c:$curb-$cure
				if {$e > $cure} {
					lappend result $c:$cure-$e
				}
				set cur [regions_insert_next $inserts pos curc curb cure]
				set region [regions_insert_next $regions regionpos c b e]
			} else {
				if {$curb > $e} {
					# puts "no overlap, not at rdna region yet, lappend region"
					lappend result $region
					set region [regions_insert_next $regions regionpos c b e]
				} elseif {$cure < $b} {
					# puts "no overlap, but rdna region is passed"
					lappend result $cur
					set cur [regions_insert_next $inserts pos curc curb cure]
				} elseif {$e > $curb} {
					# puts "overlap"
					if {$b < $curb} {
						lappend result $c:$b-$curb
					}
					lappend result $c:$curb-$cure
					while {$cure > $e} {
						set region [regions_insert_next $regions regionpos c b e]
					}
					if {$b < $cure} {
						if {$cure < $e} {
							set b $cure
							set region $c:$b-$e
						} else {
							set region [regions_insert_next $regions regionpos c b e]
						}
					}
					set cur [regions_insert_next $inserts pos curc curb cure]
				}
			}
		} elseif {$region eq "" || ($cur ne "" && [bsort [list $cur $region]] eq [list $cur $region])} {
			lappend result $cur
			set cur [regions_insert_next $inserts pos curc curb cure]
		} else {
			lappend result $region
			set region [regions_insert_next $regions regionpos c b e]
		}
	}
	return $result
}

proc regions_skip {region skipregions} {
	if {$region in $skipregions} {return 1}
	foreach skipregion $skipregions {
		if {[regexp ^$skipregion\[:_\ -\] $region]} {return 1}
	}
	return 0
}

proc getorganelles {refseq organelles} {
	if {![llength $organelles]} {
		global cache_organelles
		if {![info exists cache_organelles)]} {
			set cache_organelles {}
			set ofile [gzfile [refdir $refseq]/extra/reg_*_organelles.tsv]
			if {[file exists $ofile]} {
				foreach o [read_tsv $ofile chromosome] {
					lappend cache_organelles $o
				}
			}
		}
		return $cache_organelles
	}
	return $organelles
}

proc getrDNA {refseq rDNA} {
	if {![llength $rDNA]} {
		global cache_rDNA
		if {![info exists cache_rDNA]} {
			set cache_rDNA {}
			set ofile [gzfile [refdir $refseq]/extra/reg_*_rDNA.tsv]
			if {[file exists $ofile]} {
				set f [gzopen $ofile]
				set header [tsv_open $f]
				if {$header ne {chromosome begin end}} {error "$ofile does not have fields: chromosome begin end"}
				while {[gets $f line] != -1} {
					foreach {c b e} [split $line \t] break
					incr b
					lappend cache_rDNA $c:$b-$e
				}
				gzclose $f
			}
		}
		return $cache_rDNA
	}
	return $rDNA
}

proc regions_organelle {refseq organelles region} {
	set organelles [getorganelles $refseq $organelles]
	if {$region in $organelles} {return 1}
	foreach o $organelles {
		if {[regexp ^$o\[:_\ -\] $region]} {return 1}
	}
	return 0
}

proc regions_rDNA {refseq rDNA region} {
	set rDNA [getrDNA $refseq $rDNA]
	if {$region in $rDNA} {return 1}
	foreach {c b e} {{} {} {}} break
	foreach {c b e} [split $region -:] break
	foreach line $rDNA {
		foreach {rc rb re} {{} {} {}} break
		foreach {rc rb re} [split $line -:] break
		if {$c ne $rc} continue
		if {$rb eq "" && $re eq ""} {return 1}
		if {$b >= $rb && $e <= $re} {return 1}
	}
	return 0
}
