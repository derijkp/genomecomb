package require BioTcl

proc liftover_correctline {line ref altpos alts sposs slist} {
	if {$altpos != -1} {lset line $altpos [join $alts ,]}
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
	set regionfile {}
	set correctvariants 1
	unset -nocomplain split
	cg_options liftover args {
		-regionfile - -r {
			set regionfile $value
		}
		-correctvariants - -c {
			set correctvariants $value
		}
		-split - -s {
			set split $value
		}
	} {varfile resultfile liftoverfile} 3
	if {[file exists $resultfile]} {
		error "file $resultfile already exists, format is (now): cg liftover varfile resultfile liftoverfile"
	}
	if {[file isdir $varfile]} {
		cg_liftoversample $varfile $resultfile $liftoverfile {*}$args
	}
	set unmappedfile $resultfile.unmapped
	#
	catch {gzclose $f} ; catch {gzclose $fl} ; catch {gzclose $fc} ; catch {gzclose $freg} ; catch {close $o} ; catch {close $ou}
	#
	# open liftover file
	set fl [openliftoverfile $liftoverfile lheader oldrefname newrefname]
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
	set cinfo [comment2dict $comment]
	if {![info exists split]} {
		if {[dict exists $cinfo split]} {
			set split [dict get $cinfo split]
		} else {
			set split 1
		}
	} else {
		if {[dict exists $cinfo split]} {
			if {[dict get $cinfo split] ne $split} {
				exiterror "option -split $split was given, but file $varfile is indicated as [dict get $cinfo split]"
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
	#
	# open regionfile
	if {$regionfile ne ""} {
		set doregions 1
		set regtemplate [list_fill [llength $newheader] ?]
		lset regtemplate $typepos snp
		set freg [gzopen $regionfile]
		set regheader [tsv_open $freg]
		set regposs [tsv_basicfields $header]
		set curreg [get_region $freg $regposs]
		set samples [samples $header]
		set regsampleposs {}
		foreach sample $samples {
			set pos [lsearch $regheader $sample]
			if {$pos == -1} {
				set pos [lsearch $regheader sreg-$sample]
				if {$pos == -1} {
					set pos [lsearch $regheader reg-$sample]
					if {$pos == -1} {error "regionfile $regionfile does not contain region information for sample $sample in varfile"}
				}
			}
			lappend regsampleposs $pos
		}
	} else {
		set doregions 0
	}
	# go
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
	}
	if {[llength $aposs] && $altpos == -1} {
		set doalt 1
		set altpos [llength $header]
		lappend header alt
		lset fposs 5 $altpos
		lappend newheader alt
	} else {
		set doalt 0
	}
	lappend newheader ${oldrefname}_chromosome ${oldrefname}_begin ${oldrefname}_end ${oldrefname}_ref
	#
	# open resultfile
	set o [open $resultfile.temp w]
	dict set cinfo liftover_source $varfile
	dict set cinfo liftover $liftoverfile
	if {[dict exists $cinfo ref]} {
		set oldref [dict get $cinfo ref]
	} else {
		set oldref $oldrefname
	}
	dict set cinfo split $split
	dict set cinfo oldref $oldref
	dict set cinfo ref $newrefname
	puts -nonewline $o [dict2comment $cinfo]
	puts $o [join $newheader \t]
	# open unmappedfile
	set ou [open $unmappedfile.temp w]
	puts -nonewline $ou $comment
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
		if {$doalt} {
			set as [list_sub $line $aposs]
			set alt [join [list_remove [list_remdup $as] $ref] ,]
			lset floc 5 $alt
			lappend line $alt
		}
		set loc [lrange $floc 0 2]
		lappend line ${chromosome} ${begin} ${end} $ref
		# if {$strand ne ""} {append before -$strand}
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
				if {$doregions && $ccomp != 0} {
					set regcomp 0
					while 1 {
						set regcomp [reg_compare $cloc $curreg]
						if {$regcomp <= 0} {
							# variant is in or before region
							break
						}
						if {[eof $freg]} break
						set curreg [get_region $freg $regposs]
					}
					if {$regcomp == 0} {
						foreach {cchromosome cbegin cend csrcref cdestchromosome cdestbegin cdestend destref destcomplement} $cline break
						set temp $regtemplate
						lset temp $chrpos $cdestchromosome
						lset temp $beginpos $cdestbegin
						lset temp $endpos $cdestend
						lset temp $refpos $destref
						if {$altpos != -1} {lset temp $altpos $csrcref}
						lappend temp $cchromosome $cbegin $cend $csrcref
						foreach {a1pos a2pos seqpos zygpos} $sposs sreg [list_sub $curreg $regsampleposs] {
							if {$sreg} {
								lset temp $a1pos $csrcref
								lset temp $a2pos $csrcref
								lset temp $seqpos v
								lset temp $zygpos m
							} else {
								lset temp $seqpos u
								lset temp $zygpos u
							}
						}
						puts $o [join $temp \t]
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
					if {$refpos != -1} {lset line $refpos $ref}
					if {$altpos != -1} {lset line $altpos [join $alts ,]}
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
					foreach {cchromosome cbegin cend csrcref cdestchromosome cdestbegin cdestend cdestref cdestcomplement} $cline break
					set spos [expr {$cbegin-$begin}]
					set newref [string replace $newref $spos $spos $cdestref]
				}
				if {$newref ne $ref} {
					if {$refpos != -1} {lset line $refpos $newref}
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

	catch {gzclose $f} ; catch {closeliftoverfile $fl} ; catch {gzclose $fc} ; catch {close $o} ; catch {close $ou}
	#
	# sort result
	set sortfields [list_sub $header [list_remove [lrange $fposs 0 5] -1]]
	lappend sortfields ${oldrefname}_chromosome ${oldrefname}_begin ${oldrefname}_end ${oldrefname}_ref
	set fields [list_union [list_sub $newheader [list_remove $fposs -1]] $newheader]
	cg select -f $fields -s $sortfields $resultfile.temp $resultfile.temp2
	file rename -force $resultfile.temp2 $resultfile
	file delete -force $resultfile.temp
	#
	# rename result, cleanup
	#
	file rename -force $unmappedfile.temp $unmappedfile
}
