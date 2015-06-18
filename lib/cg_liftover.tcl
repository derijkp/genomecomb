package require BioTcl

proc liftover_correctline {line ref altpos alts sposs slist} {
	lset line $altpos [join $alts ,]
	foreach {a1pos a2pos seqpos zygpos} $sposs {a1 a2 seq zyg} $slist {
		if {$a1pos == -1 || $a2pos == -1} continue
		if {$seq eq "u"} {
			set zyg u
		} else {
			set zyg [zyg $a1 $a2 $ref $alts]
			if {$zyg in "r o"} {set seq r} else {set seq v}
		}
		if {$zygpos != -1} {
			lset line $zygpos $zyg
		}
		if {$seqpos != -1} {
			lset line $seqpos $seq
		}
	}
	return $line
}

proc cg_liftover {args} {
	set pos 0
	set dbdir {}
	set regionfile {}
	set correctvariants 1
	unset -nocomplain split
	foreach {key value} $args {
		switch -- $key {
			-dbdir {
				set dbdir $value
			}
			-regionfile - -r {
				set regionfile $value
			}
			-correctvariants - -c {
				set correctvariants $value
			}
			-split - -s {
				set split $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 3)} {
		errorformat liftover
		exit 1
	}
	foreach {varfile resultfile liftoverfile} $args break
	if {[file exists $resultfile]} {
		error "file $resultfile already exists, format is (now): cg liftover varfile resultfile liftoverfile"
	}

	set liftoverfile [liftoverfile $liftoverfile]
	if {[file isdir $varfile]} {
		cg_liftoversample {*}$args
	}
	set unmappedfile $resultfile.unmapped
	#
	catch {gzclose $f} ; catch {gzclose $fl} ; catch {gzclose $fc} ; catch {close $o} ; catch {close $ou}
	#
	# open liftover file
	set fl [gzopen $liftoverfile]
	set lheader [tsv_open $fl comment]
	set lposs [list_cor $lheader {chromosome begin end strand destchromosome destbegin destend deststrand}]
	set lline [list_sub [split [gets $fl] \t] $lposs]
	set fromloc [lrange $lline 0 2]
	foreach {srcchromosome srcbegin srcend srcstrand destchromosome destbegin destend deststrand} $lline break
	if {-1 in $lposs} {exiterror "error in liftoverfile ($liftoverfile): wrong header"}
	#
	# open refchanges
	if {$correctvariants} {
		set refchangesfile [file root $liftoverfile].refchanges.tsv
		set refchangesfile [gzfiles $refchangesfile]
		if {$refchangesfile eq ""} {
			exiterror "option correctvariants is 1, but could not find the needed refchanges file $refchangesfile\nThis can be generated using the command:\ncg liftfindchanges srcgenome.ifas destgenome.ifas $liftoverfile > $refchangesfile"
		}
		set fc [gzopen $refchangesfile]
		set cheader [tsv_open $fc]
		if {$cheader ne "chromosome begin end ref destchromosome destbegin destend destref destcomplement"} {
			exiterror "incorrect header for refchangesfile $refchangesfile"
		}
		# clist will gather potential overlapping refchanges
		set clist {}
	}
	#
	# open varfile
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set commentd [comment2dict $comment]
	if {![info exists split]} {
		if {[dict exists $commentd split]} {
			set split [dict get $commentd split]
		} else {
			set split 1
		}
	} else {
		if {[dict exists $commentd split]} {
			if {[dict get $commentd split] ne $split} {
				exiterror "option -split $split was given, but file $varfile is indicated as [dict get $commentd split]"
			}
		}
	}
	# next line here to give error if first 3 are not present
	set poss [tsv_basicfields $header 3]
	set fposs [tsv_basicfields $header 6 0]
	lappend fposs [lsearch $header strand]
	foreach {chrpos beginpos endpos typepos refpos altpos strandpos} $fposs break
	set strand {}
	set newheader $header
	lappend newheader beforeliftover
	set aposs {}
	set sposs {}
	if {$correctvariants} {
		set samples [samples $header]
		if {[llength $samples]} {
			foreach sample $samples {
				lappend sposs {*}[list_cor $header [list alleleSeq1-$sample alleleSeq2-$sample sequenced-$sample zyg-$sample]]
				lappend aposs {*}[lrange $sposs end-3 end-2]
			}
		} else {
			set tempposs [list_cor $header [list alleleSeq1 alleleSeq2 sequenced zyg]]
			if {[llength [list_remove $tempposs -1]]} {
				lappend sposs {*}$tempposs
				lappend aposs {*}[lrange $sposs end-3 end-2]
			}
		}
		if {[llength $aposs]} {
			set doalt 1
		} else {
			set doalt 0
		}
	}
	#
	# open resultfile
	set o [open $resultfile.temp w]
	puts -nonewline $o $comment
	puts $o "#liftover_source\t$varfile"
	puts $o "#liftover\t$liftoverfile"
	puts $o [join $newheader \t]
	# open unmappedfile
	set ou [open $unmappedfile.temp w]
	puts $ou "#liftover_source\t$varfile"
	puts $ou "#liftover_unmapped\t$liftoverfile"
	puts $ou [join $header \t]
	set ldone 0
	set todo {}
	set doneloc {}
	set donetype {}

	while 1 {
		if {[gets $f oline] == -1} break
		set line [split $oline \t]
		set floc [list_sub $line $fposs]
		foreach {chromosome begin end type ref alt strand} $floc break
		set loc [lrange $floc 0 2]
		set before ${chromosome}-${begin}-${end}
		if {$strand ne ""} {append before -$strand}
		lappend line $before
		if {$correctvariants} {
			# adapt clist to current var
			set newclist {}
			# coverlaps will gather refchanges overlapping with current var
			set coverlaps {}
			set ccomp 0
			# check if refchanges in clist can be dropped (come before current var pos)
			foreach cline $clist {
				set cloc [lrange $cline 0 2]
				set ccomp [reg_compare $cloc $loc]
				if {$ccomp >= 0} {
					lappend newclist $cline
					if {$ccomp == 0} {
						lappend coverlaps $cline
					}
				}
			}
			# check if new refchanges should be added to clist
			# stops if new refchange as after current var pos
			while {$ccomp <= 0 && ![eof $fc]} {
				if {[gets $fc cline] == -1} break
				set cline [split $cline \t]
				set cloc [lrange $cline 0 2]
				set ccomp [reg_compare $cloc $loc]
				if {$ccomp >= 0} {
					lappend newclist $cline
					if {$ccomp == 0} {
						lappend coverlaps $cline
					}
				}
			}
			set clist $newclist
		}
		# find transfer region applicable to current var
		while 1 {
			set comp [reg_compare $fromloc $loc]
			if {$comp >= 0} {
				break
			}
			if {$ldone} break
			if {[gets $fl lline] == -1} {
				set ldone 1
				break
			}
			set lline [list_sub [split $lline \t] $lposs]
			foreach {srcchromosome srcbegin srcend srcstrand destchromosome destbegin destend deststrand} $lline break
			set fromloc [lrange $lline 0 2]
		}
		if {$comp == 0 && $begin >= $srcbegin && $end <= $srcend} {
			# applicable transfer region found -> remap var
			set alts [split $alt ,]
			set slist [list_sub $line $sposs]
			set alist [list_sub $line $aposs]
			if {$srcstrand eq $deststrand} {
				set ustrand $strand
				set ubegin [expr {$begin + $destbegin - $srcbegin}]
				set uend [expr {$end + $destbegin - $srcbegin}]
			} else {
				if {$strand eq "+"} {set ustrand "-"} else {set ustrand "+"}
				set uend [expr {$destend - $begin + $srcbegin}]
				set ubegin [expr {$destend - $end + $srcbegin}]
				set ref [seq_complement $ref]
				set alts [lsort -dict [seq_complement $alts]]
				if {$correctvariants} {
					lset line $refpos $ref
					lset line $altpos [join $alts ,]
					foreach {a1pos a2pos seqpos zygpos} $sposs {a1 a2 seq zyg} $slist {
						set a1 [seq_complement [lindex $line $a1pos]]
						set a2 [seq_complement [lindex $line $a2pos]]
						lset line $a1pos $a1
						lset line $a2pos $a2
					}
				}
				set slist [list_sub $line $sposs]
				set alist [list_sub $line $aposs]
			}
			lset line $chrpos $destchromosome
			lset line $beginpos $ubegin
			lset line $endpos $uend
			if {$strandpos != -1} {
				lset line $strandpos $strand
			}
			if {$correctvariants && [llength $coverlaps]} {
				set newref $ref
				foreach cline $coverlaps {
					foreach {cchromosome cbegin cend cstrand srcref temp temp destref destcomplement} $cline break
					set spos [expr {$cbegin-$begin}]
					set newref [string replace $newref $spos $spos $destref]
				}
				if {$newref ne $ref} {
					lset line $refpos $newref
					if {!$split} {
						set alts [list_remove $alts $newref]
						lappend alts $ref
						set alts [lsort -dict $alts]
						puts $o [join [liftover_correctline $line $newref $altpos $alts $sposs $slist] \t]
					} else {
						if {[llength $alts] > 1} {error "file has multiple alternative alleles (but should be split)"}
						set alt [lindex $alts 0]
						# alts should be only 1 long
						if {(![llength $alist] || $ref in $alist) && $type ne $donetype || $doneloc ne $loc} {
							# output (but only once) line with old ref as alt allele
							set doneloc $loc
							set donetype $type
							puts $o [join [liftover_correctline $line $newref $altpos $ref $sposs $slist] \t]
						}
						if {$alt ne $newref} {
							puts $o [join [liftover_correctline $line $newref $altpos $alt $sposs $slist] \t]
						}
					}
					set ref $newref
				} else {
					puts $o [join $line \t]
				}
			} else {
				puts $o [join $line \t]
			}
		} else {
			# no applicable transfer region found -> unmapped
			puts $ou $oline
		}
	}
	catch {gzclose $f} ; catch {gzclose $fl} ; catch {gzclose $fc} ; catch {close $o} ; catch {close $ou}
	#
	# sort result
	set sortfields [list_sub $header [lrange $fposs 0 5]]
	lappend sortfields beforeliftover
	cg select -s {chromosome begin end type ref alt beforeliftover} $resultfile.temp $resultfile.temp2
	file rename -force $resultfile.temp2 $resultfile
	file delete -force $resultfile.temp
	#
	# rename result, cleanup
	#
	file rename -force $unmappedfile.temp $unmappedfile
}
