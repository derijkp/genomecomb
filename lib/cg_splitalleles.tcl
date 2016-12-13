proc cg_splitalleles {args} {
#	set pos 0
#	while 1 {
#		set key [lindex $args $pos]
#		switch -- $key {
#			-o {
#				incr pos
#				set outdir [file_absolute [lindex $args $pos]]
#				incr pos
#			}
#			-- break
#			default {
#				break
#			}
#		}
#	}
#	set args [lrange $args $pos end]
	foreach {file outfile} {{} {}} break
	cg_options splitalleles args {
	} {file outfile} 0 2
	if {$file eq ""} {
		set f stdin
	} else {
		set f [gzopen $file]
	}
	if {$outfile eq ""} {
		set o stdout
	} else {
		set o [open $outfile w]
	}
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 6 0]
	set apos [lindex $poss 5]
	if {$apos == -1} {error "alt filed not found"}
	set rpos [lindex $poss 4]
	set samples [samples $header]
	set sposs {}
	foreach sample $samples {
		lappend sposs {*}[list_cor $header [list alleleSeq1-$sample alleleSeq2-$sample sequenced-$sample zyg-$sample]]
	}
	puts -nonewline $o $comment
	puts $o [join $header \t]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set ref [lindex $line $rpos]
		set alt [lindex $line $apos]
		set alleles [split $alt ,]
		# redundancy in alleles is an error, but I have come across it, so dedup
		set alleles [list_remdup $alleles]
		lset line $apos [join $alleles ,]
		if {![llength $alleles]} {set alleles [list {}]}
		unset -nocomplain a
		if {[llength $alleles] == 1} {
			set a([lindex $alleles 0]) $line
		} else {
			set alen [llength $alleles]
			foreach vs $line {
				set svs [split $vs ,]
				if {[llength $svs] <= 1 || [llength $svs] != $alen} {
					foreach allele $alleles {
						lappend a($allele) $vs
					}
				} else {
					foreach allele $alleles v $svs {
						lappend a($allele) $v
					}
				}
			}
			# foreach all $alleles {puts $o "$all [lrange $a($all) 0 8]"}
		}
		foreach alt [ssort -natural $alleles] {
			set line $a($alt)
			foreach {a1pos a2pos seqpos zygpos} $sposs {
				set a1 [lindex $line $a1pos]
				set a2 [lindex $line $a2pos]
				set seq [lindex $line $seqpos]
				if {$seq eq "u"} {
					set zyg u
				} else {
					set zyg [zyg $a1 $a2 $ref $alt]
					if {$zyg in "r o"} {set seq r} else {set seq v}
				}
				if {$zygpos != -1} {
					lset line $zygpos $zyg
				}
				if {$seqpos != -1} {
					lset line $seqpos $seq
				}
			}
			set a($alt) $line
		}
		foreach allele [ssort -natural $alleles] {
			puts $o [join $a($allele) \t]
		}
	}
	gzclose $f
}
