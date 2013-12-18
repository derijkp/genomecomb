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
	if {[llength $args] != 1} {
		error "format is: cg splitalleles file"
		exit 1
	}
	foreach file $args break
	set f [gzopen $file]
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header]
	set apos [lindex $poss 5]
	set rpos [lindex $poss 4]
	set samples [samples $header]
	set sposs {}
	foreach sample $samples {
		lappend sposs {*}[list_cor $header [list alleleSeq1-$sample alleleSeq2-$sample sequenced-$sample]]
	}
	if {[string length $comment]} {
		puts $comment
	}
	puts [join $header \t]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set ref [lindex $line $rpos]
		set alleles [split [lindex $line $apos] ,]
		if {![llength $alleles]} {set alleles [list {}]}
		if {[llength $alleles] == 1} {
			set a([lindex $alleles 0]) $line
		} else {
			unset -nocomplain a
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
			# foreach all $alleles {puts "$all [lrange $a($all) 0 8]"}
		}
		foreach alt [ssort -natural $alleles] {
			set line $a($alt)
			foreach {a1pos a2pos seqpos} $sposs {
				set a1 [lindex $line $a1pos]
				set a2 [lindex $line $a2pos]
				set seq [lindex $line $seqpos]
				if {$seq eq "u"} {
					set seq u
				} elseif {$a1 eq $alt} {
					if {$a2 eq $alt} {
						set seq m
					} elseif {$a2 eq $ref} {
						set seq t
					} else {
						set seq c
					}
				} elseif {$a2 eq $alt} {
					if {$a2 eq $ref} {
						set seq t
					} else {
						set seq c
					}
				} elseif {$a1 eq $ref && $a2 eq $ref} {
					set seq r
				} else {
					set seq o
				}
				if {$seqpos != -1} {
					lset line $seqpos $seq
				}
			} 
		} 
		foreach allele [ssort -natural $alleles] {
			puts [join $a($allele) \t]
		}
	}
	close $f
}
