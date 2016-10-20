#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multicompar_reannot_find {basedir sample args} {
	if {![llength $args]} {set args [list {}]}
	set sampledir [lindex [split $sample -] end]
	if {![llength $args]} {
		set args [list {}]
		set samplelist [list $sample $sampledir]
	} else {
		set samplelist [list $sample $sampledir {}]
	}
	foreach startdir [list $basedir/samples $basedir [file dir $basedir]/samples [file dir $basedir]] {
		foreach usesample $samplelist {
			foreach pattern $args {
				set test [gzfile [file join $startdir $usesample $pattern]]
				if {[file exists $test]} {
					return $test
				}
			}
		}
	}
	return {}
}

proc multicompar_reannot_getline {f sampleaVar poss} {
	global reannot_cache reannot_prev
	upvar $sampleaVar samplea
	if {[llength $reannot_cache] > 1} {
		return [list_shift reannot_cache]
	}
	while {1} {
		set line [split [gets $f] \t]
		if {[llength $line]} {
			set cur [list_sub $line $poss]
			set d [lloc_compare $reannot_prev $cur]
			if {$d < 0} {
				break
				set reannot_prev $cur
			} elseif {$d == 0} {
				lappend reannot_cache $line
			} else {
				error "Cannot reannot because the file is not correctly sorted (sort correctly using \"cg select -s -\")"
			}
		} elseif {[eof $f]} break
	}
	if {[llength $reannot_cache] == 1} {
		set result [lindex $reannot_cache 0]
		set reannot_cache [list $line]
		set reannot_prev [list_sub $line $poss]
		return $result
	} else {
		foreach pos $samplea(aposs) {
			set list [list_subindex $reannot_cache $pos]
			set slist [list_remove $list ?]
			if {![llength $slist]} continue
			set a [lindex $slist 0]
			set i 0
			foreach ta $list {
				if {$ta eq "?"} {
					lset reannot_cache $i $pos $a
				}
				incr i
			}
		}
		lappend reannot_cache $line
		set reannot_prev [list_sub $line $poss]
		return [list_shift reannot_cache]
	}
}

proc multicompar_reannot {compar_file {force 0} {regonly 0} {skipincomplete 0} {range {}} {sampleaVar {}}} {
# putsvars compar_file force regonly skipincomplete range sampleaVar
	if {$sampleaVar ne ""} {upvar $sampleaVar samplea}
	set compar_file [file_absolute $compar_file]
	set basedir [file dir $compar_file]
	catch {close $f}; catch {close $o}
	set f [gzopen $compar_file]
	set header [tsv_open $f]
	set pos -1
	# allow reuse of cached samplea for paging
	# unset -nocomplain samplea
	if {![info exists samplea(samples)]} {
		set samples {}
		foreach field $header {
			incr pos
			set splitpos [string first - $field]
			if {$splitpos != -1} {
				set temp [string range $field 0 [expr {$splitpos-1}]]
				incr splitpos
				set sample [string range $field $splitpos end]
				lappend samplea(poss,$sample) $pos
				lappend samplea(fields,$sample) $temp
				lappend samples $sample
			}
		}
		set samplea(samples) [list_remdup $samples]
	}
	set samples $samplea(samples)
	if {[llength $range]} {
		set samples [lrange $samples {*}$range]
		if {![llength $samples]} {return {}}
	}
	set referencepos [lsearch $header reference]
	if {$referencepos == -1} {set referencepos [lsearch $header ref]}
	putslog "file contains samples: $samples"
	set samplea(aposs) {}
	foreach sample $samples {
		set poss [tsv_basicfields $header 6 0]
		set samplea(ref,$sample) [lindex $poss 4]
		set samplea(alt,$sample) [lindex $poss 5]
		set samplea(a1,$sample) [lsearch $header alleleSeq1-$sample]
		lappend samplea(aposs) $samplea(a1,$sample)
		set samplea(a2,$sample) [lsearch $header alleleSeq2-$sample]
		lappend samplea(aposs) $samplea(a2,$sample)
		set samplea(rpos,$sample) [lsearch $header refscore-$sample]
		set samplea(cpos,$sample) [lsearch $header coverage-$sample]
		set samplea(seq,$sample) [lsearch $header sequenced-$sample]
		set samplea(zyg,$sample) [lsearch $header zyg-$sample]
		set samplea(dir,$sample) [multicompar_reannot_find $basedir $sample]
		set samplea(regionfile,$sample) [multicompar_reannot_find $basedir $sample sreg-$sample.tsv]
		set samplea(annot_regionfile,$sample) [file exists $samplea(regionfile,$sample)]
		set samplea(annot_coverage,$sample) 0
		set okinfo 0
		if {[file exists $samplea(dir,$sample)/allpos]} {
			set samplea(type,$sample) rtg
			annot_rtg_init $samplea(dir,$sample)
			set samplea(rtgposs,$sample) {}
			foreach field {
				alleleSeq1 alleleSeq2 posterior coverage correction
				numA numC numG numT percA percC percG percT nonidentityposterior
			} {
				lappend samplea(rtgposs,$sample) [lsearch $header ${field}-$sample]
			}
			set okinfo 1
		} elseif {$samplea(annot_regionfile,$sample)} {
			set samplea(type,$sample) cg
			if {!$regonly} {
				putslog "Using $samplea(dir,$sample) to reannot"
				annot_coverage_init $samplea(dir,$sample) $sample
				set samplea(annot_coverage,$sample) 1
			}
			putslog "Using $samplea(regionfile,$sample) to reannot"
			annot_region_init $samplea(regionfile,$sample)
			set okinfo 1
		}
		set samplea(varall,$sample) [multicompar_reannot_find $basedir $sample varall-$sample.tsv]
		if {!$regonly && $samplea(varall,$sample) ne ""} {
			set samplea(type,$sample) varall
			putslog "Using $samplea(varall,$sample) to reannot"
			annot_varall_init $samplea(varall,$sample) $sample $header
			set okinfo 1
		}
		if {!$okinfo} {
			if {!$skipincomplete} {
				error "no sorted region file (sreg-$sample.tsv) or allpos dir (for rtg) found: not properly processed sample"
			} else {
				putslog "warning: no sorted region file (sreg-$sample.tsv) or allpos dir (for rtg) found: not properly processed sample"
			}
			set samples [list_remove $samples $sample]
		}
		set samplea(todo,$sample) {}
		foreach temp {refcons nocall cluster} {
			set regfile [multicompar_reannot_find $basedir $sample reg_${temp}-$sample.tsv]
			if {$regfile ne "" && [inlist $samplea(fields,$sample) $temp]} {
				lappend samplea(todo,$sample) [list [lsearch $header ${temp}-$sample] 1 $regfile]
			}
		}
		list_foreach {field value regfile} $samplea(todo,$sample) {
			annot_region_init $regfile
		}
	}
	# start processing
	set o [open $compar_file.temp w]
	puts $o [join $header \t]
	set poss [tsv_basicfields $header 4]
	set num 0
	set line [split [gets $f] \t]
	set ::reannot_cache [list $line]
	set ::reannot_prev [list_sub $line $poss]

	while {1} {
		incr num
		if {![expr $num%10000]} {putslog $num}
		set line [multicompar_reannot_getline $f samplea $poss]
		if {![llength $line]} break
		set reference [lindex $line $referencepos]
		set loc [list_sub $line $poss]
		foreach {chr begin end} $loc break
		foreach sample $samples {
			#if {!$force && ([lindex $line $samplea(a1,$sample)] ne "?")} continue
			list_foreach {field value regfile} $samplea(todo,$sample) {
				if {[lindex $line $field] == "-"} continue
				if {!$force && [lindex $line $field] != "?"} continue
				set r [annot_region_get $regfile $chr $begin $end]
				if {$r} {lset line $field $value} else {lset line $field {}}
			}
			if {$samplea(annot_coverage,$sample) && ($force || ([lindex $line $samplea(rpos,$sample)] eq "?") || ([lindex $line $samplea(cpos,$sample)] eq "?"))} {
				# from cg, combined refscore and coverage (would be better separate if gotten from bcols, but can also be from one tsv file)
				foreach {r c} [annot_coverage_get $samplea(dir,$sample) $sample $chr $begin] break
				if {$samplea(rpos,$sample) != -1} {lset line $samplea(rpos,$sample) $r}
				lset line $samplea(cpos,$sample) $c
			}
			if {$samplea(type,$sample) eq "varall"} {
				# varall is used, region files are also used if present
				annot_varall_annot samplea $sample $loc $force line
			} elseif {$samplea(type,$sample) eq "cg"} {
				# region files only for seq and zyg
				set seq [lindex $line $samplea(seq,$sample)]
				set zyg [lindex $line $samplea(zyg,$sample)]
				if {$seq eq "?" || $zyg eq "?"} {
					set r [annot_region_in $samplea(regionfile,$sample) $chr $begin $end]
					set a1 [lindex $line $samplea(a1,$sample)]
					set a2 [lindex $line $samplea(a2,$sample)]
					if {$r} {
						set a1 [lindex $line $samplea(a1,$sample)]
						set a2 [lindex $line $samplea(a2,$sample)]
						if {[inlist {- ?} $a1]} {
							set a1 $reference
							lset line $samplea(a1,$sample) $reference
						}
						if {[inlist {- ?} $a2]} {
							set a2 $reference
							lset line $samplea(a2,$sample) $reference
						}
						if {$samplea(alt,$sample) != -1} {
							set alt [lindex $line $samplea(alt,$sample)]
						} else {
							set alt ?
						}
						set zyg [zyg $a1 $a2 $reference $alt]
						if {$zyg in "r o"} {set seq r} else {set seq v}
						lset line $samplea(seq,$sample) $seq
						if {$samplea(zyg,$sample) != -1} {lset line $samplea(zyg,$sample) $zyg}
					} else {
						lset line $samplea(seq,$sample) u
						if {$samplea(zyg,$sample) != -1} {lset line $samplea(zyg,$sample) u}
						# if {[inlist {- ?} $a1]} {lset line $samplea(a1,$sample) -}
						# if {[inlist {- ?} $a2]} {lset line $samplea(a2,$sample) -}
					}
				}
			} elseif {$samplea(type,$sample) eq "rtg"} {
				set sub [list_sub $line $samplea(rtgposs,$sample)]
				if {!$force && ![inlist [list_sub $line $samplea(rtgposs,$sample)] ?]} continue
				set temp [annot_rtg_get $samplea(dir,$sample) $chr $begin]
				set rtgdata [lrange $temp 4 end-1]
				if {[llength $rtgdata] < [llength $samplea(rtgposs,$sample)]} {
					lset line $samplea(seq,$sample) u
					foreach pos $samplea(rtgposs,$sample) v $rtgdata {
						if {$pos == -1} continue
						lset line $pos -
					}
				} elseif {[llength $rtgdata] > [llength $samplea(rtgposs,$sample)]} {
					error "error in rtg data\nsample=$sample\nline=$line"
				} else {
					set coverage [lindex $rtgdata 3]
					if {$coverage < 10} {
						lset line $samplea(seq,$sample) u
					} else {
						lset line $samplea(seq,$sample) r
					}
					foreach pos $samplea(rtgposs,$sample) v $rtgdata {
						if {$pos == -1} continue
						lset line $pos $v
					}
				}
			}
		}
		puts $o [join $line \t]
	}

	close $o
	gzclose $f
	foreach sample $samples {
		list_foreach {field value regfile} $samplea(todo,$sample) {
			annot_region_close $regfile
		}
		if {$samplea(varall,$sample) ne ""} {
			if {!$regonly} {
				annot_varall_close $samplea(varall,$sample) $sample
			}
		}
		if {$samplea(type,$sample) eq "cg"} {
			if {!$regonly} {
				annot_coverage_close $samplea(dir,$sample) $sample
			}
			annot_region_close $samplea(regionfile,$sample)
		} elseif {$samplea(type,$sample) eq "rtg"} {
			annot_rtg_close $samplea(dir,$sample)
		}
	}
	file rename -force $compar_file $compar_file.old
	file rename -force $compar_file.temp $compar_file
	return $samples
}


proc cg_multicompar_reannot {args} {
	set force 0
	set regonly 0
	set skipincomplete 0
	set paged 0
	set pagedstart 0
	set pos 0
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-force {
				incr pos
				set force [lindex $args $pos]
			}
			-regonly {
				incr pos
				set regonly [lindex $args $pos]
			}
			-skipincomplete {
				incr pos
				set skipincomplete [lindex $args $pos]
			}
			-paged {
				incr pos
				set paged [lindex $args $pos]
			}
			-pagedstart {
				incr pos
				set pagedstart [lindex $args $pos]
			}
			-- break
			default {
				if {[string index $key 0] eq "-"} {error "unknown option \"$key\""}
				break
			}
		}
		incr pos
	}
	set args [lrange $args $pos end]
	if {[llength $args] < 1} {
		errorformat multicompar_reannot
		exit 1
	}
	set compar_file [list_shift args]
	foreach option $args {
		switch $option {
			force {set force 1}
			regonly {set regonly 1}
			skipincomplete {set skipincomplete 1}
			paged {set paged 1000}
			default {error "unrecognized option $option"}
		}
	}
	putslog "Reannotating $compar_file"
	if {!$paged} {
		multicompar_reannot $compar_file $force $regonly $skipincomplete
	} else {
		set s $pagedstart
		unset -nocomplain samplea
		while 1 {
			set e [expr {$s + $paged-1}]
			set samples [multicompar_reannot $compar_file $force $regonly $skipincomplete [list $s $e] samplea]
			if {![llength $samples]} break
			putslog "range $s-$e done"
			incr s $paged
		}
	}
}
