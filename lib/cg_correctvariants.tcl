#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

if 0 {
	cd /complgen/projects/dlb1
	set file oldcmt71_compar.tsv.old
}

proc cg_correctvariants_alts {type list gref} {
	set list [list_remdup $list]
	set alts [list_remove $list - ? N $gref {}]
	if {([llength $alts] == 0) && ($type ne "del")} {set alts ?}
	return $alts
}

proc cg_correctvariants {args} {
	set force 0
	set complement 0
	set split 1
	cg_options correctvariants args {
		-c {
			set complement $value
		}
		-f {
			set force $value
		}
		-split - -s {
			set split $value
		}
	} {file resultfile dbdir} 3 3
	set dbdir [dbdir $dbdir]
	set countdoubles 0
	set countcorrected 0
	if {[file exists $resultfile]} {error "$resultfile exists"}
	catch {close $o} ; catch	{gzclose $f}
	set f [gzopen $file]
	set header [tsv_open $f comment]
	# next line here to give error if first 3 are not present
	set fposs [tsv_basicfields $header 3]
	set fposs [tsv_basicfields $header 6 0]
	set o [open $resultfile.temp w]
	puts -nonewline $o $comment
	set nheader [list_concat {chromosome begin end type ref alt} [list_sub $header -exclude $fposs]]
	set poss [list_cor $header $nheader]
	set poss [lreplace $poss 0 5 {*}$fposs]
	set samples [samples $header]
	set aposs {}
	set sposs {}
	if {[llength $samples]} {
		foreach sample $samples {
			lappend sposs {*}[list_cor $nheader [list alleleSeq1-$sample alleleSeq2-$sample sequenced-$sample zyg-$sample]]
			lappend aposs {*}[lrange $sposs end-3 end-2]
		}
	} else {
		set tempposs [list_cor $nheader [list alleleSeq1 alleleSeq2 sequenced zyg]]
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
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	puts $o [join $nheader \t]
	set count 0
	set nextline [split [gets $f] \t]
	set nextline [list_sub $nextline $poss]
	set prevvar [lrange $nextline 0 3]
	set lines [list $nextline]
	while {![eof $f]} {
		# gather all lines on same location:
		# in case of split the may need to be resorted
		# in case of non-split the may need to be gathered
		while {![eof $f]} {
			set nextline [split [gets $f] \t]
			if {![llength $nextline]} continue
			set nextline [list_sub $nextline $poss]
			set var [lrange $nextline 0 3]
			if {$var ne $prevvar} break
			lappend lines $nextline
		}
		if {[llength $lines] > 1} {
			if {$split} {
				set lines [lsort -index 5 $lines]
			} else {
				# better would be choose best/real merge, but not now
				set line [lindex $lines 0]
				set alts {}
				foreach as [list_subindex $lines 5] {
					lappend alts {*}[split $as ,]
				}
				lset line 5 [join [bsort [list_remdup $alts]] ,]
				set lines [list $line]
			}
		}
		if {![llength $lines]} continue
		incr count
		if {$count > 1000000} {
			putslog $chr:$start-$end
			set count 0
		}
		foreach {chr start end type ref alt} [lindex $lines 0] break
		set size [expr {$end-$start}]
		if {![isint $ref] || $size <= 10} {
			if {![catch {genome_get $fg $chr $start $end} gref]} {
				set gref [string toupper $gref]
			} else {
				putslog "WARNING: Could not get reference sequence for $chr:$start-$end: not checked"
				set gref $ref
			}
		} else {
			set gref $size
		}
		if {$force == 3 || ($gref ne $ref && !($ref eq "" && $size ne ""))} {
			set resultlines {}
			set prevalt ___
			foreach line $lines {
				foreach {chr start end type ref alt} $line break
				if {$alt eq $prevalt} {incr countdoubles ; continue}
				incr countcorrected
				set prevalt $alt
				set alts [split $alt ,]
				if {$split && [llength $alts] > 1} {
					error "error: the split option is used and file contains multiallelic variants"
				}
				if {$complement && ([seq_complement $gref] eq $ref)} {
					lset line 4 $gref
					if {!$split && $doalt} {
						set altlist {}
						foreach a [list_remdup [list_sub $line $aposs]] {
							lappend altlist [seq_complement $a]
						}
						set alt [join [cg_correctvariants_alts $type $altlist $gref] ,]
					} else {
						set nalts {}
						foreach a $alts {lappend nalts [seq_complement $a]}
						set alt [join $nalts ,]
					}
					lset line 5 $alt
					foreach {a1pos a2pos seqpos zygpos} $sposs {
						if {$a1pos == -1 || $a2pos == -1} continue
						set a1 [seq_complement [lindex $line $a1pos]]
						set a2 [seq_complement [lindex $line $a2pos]]
						lset line $a1pos $a1
						lset line $a2pos $a2
						set altsa($a1) 1 ; set altsa($a2) 1
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
				} elseif {!$force} {
					error "different ref ($ref) for line (ref should be $gref):\n$line"
				} else {
					if {$force == 2} {
						putslog "different ref ($ref) for line (ref should be $gref):\n$line"
					}
					lset line 4 $gref
					if {!$split && $doalt} {
						set altlist [list_remdup [list_sub $line $aposs]]
						set alt [join [cg_correctvariants_alts $type $altlist $gref] ,]
						lset line 5 $alt
					} else {
						if {[inlist $alts $gref]} {
							set nalts [list_remove $alts $gref]
							lappend nalts $ref
							set alt [join $nalts ,]
							lset line 5 $alt
						}
					}
					foreach {a1pos a2pos seqpos zygpos} $sposs {
						if {$a1pos == -1 || $a2pos == -1} continue
						set a1 [lindex $line $a1pos]
						set a2 [lindex $line $a2pos]
						set seq [lindex $line $seqpos]
						if {$seq eq "u"} {
							set zyg u
						} else {
							set zyg [zyg $a1 $a2 $gref $alt]
							if {$zyg in "r o"} {set seq r} else {set seq v}
						}
						if {$zygpos != -1} {
							lset line $zygpos $zyg
						}
						if {$seqpos != -1} {
							lset line $seqpos $seq
						}
					}
				}
				lappend resultlines $line
			}
			if {[llength $resultlines] > 1} {
				set resultlines [lsort -index 5 $resultlines]
			}
			foreach line $resultlines {
				puts $o [join $line \t]
			}
		} else {
			set prevalt ___
			foreach line $lines {
				set alt [lindex $line 5]
				if {$alt eq $prevalt} {incr countdoubles ; continue}
				set prevalt $alt
				lset line 4 $gref
				if {!$split && $doalt} {
					lset line 5 [join [cg_correctvariants_alts $type [list_sub $line $aposs] $gref] ,]
				}
				puts $o [join $line \t]
			}
		}
		set prevvar [lrange $nextline 0 3]
		set lines [list $nextline]
	}
	close $o
	gzclose $f
	putslog "removed $countdoubles doubles"
	putslog "corrected $countcorrected refs"
	file rename -force -- $resultfile.temp $resultfile
}

proc cg_updatevarfile {args} {
	cg_correctvariants {*}$args
}
