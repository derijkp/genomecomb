#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require BioTcl

proc annotatemir_one {loc geneobj} {
	foreach {chrom snpstart snpend snptype ref alt} $loc break
	unset -nocomplain adata
	array set adata $geneobj
	set complement $adata(complement)
	set dbstart [lindex $adata(annotlist) 0 2]
	set dbend [lindex $adata(annotlist) end 3]
	set hstart $adata(start)
	set hend $adata(end)
	array set da {}
	if {$snpstart <= $hstart && $snpend >= $hend} {
		set impact GENEDEL
	} else {
		set impact {}
		set dist {}
		list_foreach {annot distback b e} $adata(annotlist) {
			if {$snpend > $b && $snpstart < $e} {
				set temp $annot
				if {$distback == 1} {
					if {!$complement} {append temp -} else {append temp +}
					append temp [expr {$e - $snpend + 1}]
				} elseif {$distback == 2} {
					if {!$complement} {append temp +} else {append temp -}
					append temp [expr {$snpstart - $b + 1}]
				} elseif {!$complement} {
					append temp +[expr {$snpstart - $b + 1}]
				} else {
					append temp +[expr {$e - $snpend + 1}]
				}
				lappend impact $temp
			}
		}
	}
	return [list [join $impact &] $adata(genename)]
}

proc annotatemir_dbopen {dbfile genecol transcriptcol} {
	set df [open_mirfile $dbfile header dposs $genecol $transcriptcol]
	# insert empty space for genobj
	set dposs [linsert $dposs 3 -1]
	while {![eof $df]} {
		set fdbline [split [gets $df] \t]
		set fdbline [list_sub $fdbline $dposs]
		foreach {dbchr dbstart dbend} $fdbline break
		if {![isint $dbstart] || ![isint $dbend]} continue
		break
	}
	set dbobj [dict create df $df dposs $dposs fdbline $fdbline prevdbloc [lrange $fdbline 0 2]]
	return $dbobj
}

proc annotatemir_dbgetlist {dbobjVar dblistVar chr start end} {
	upvar $dbobjVar dbobj
	upvar $dblistVar dblist
	set df [dict get $dbobj df]
	set prevdbloc [dict get $dbobj prevdbloc]
	set fdbline [dict get $dbobj fdbline]
	if {![llength $fdbline]} {return $dblist}
	set dposs [dict get $dbobj dposs]
	foreach {dbchr dbstart dbend} $fdbline break
	incr dbstart -2000
	incr dbend 2000
	while {![eof $df]} {
		set chrcompar [chr_compare $dbchr $chr]
		if {$chrcompar > 0} break
		if {$chrcompar == 0} {
			if {$dbstart >= $end} break
			if {$dbend >= $start} {
				lappend dblist $fdbline
				set fdbline {}
			}
		} else {
			set dblist {}
		}
		set ok 0
		while {![eof $df]} {
			set fdbline [split [gets $df] \t]
			set fdbline [list_sub $fdbline $dposs]
			foreach {dbchr dbstart dbend} $fdbline break
			if {[isint $dbstart] && [isint $dbend]} {
				set ok 1
				break
			}
		}
		#break if eof
		if {[eof $df]} break
		set dbloc [lrange $fdbline 0 2]
		if {[lloc_compare $prevdbloc $dbloc] > 0} {
			error "Cannot annotate because the database file ($dbfile) is not correctly sorted (sort correctly using \"cg select -s -\")"
		}
		set prevdbloc $dbloc
		if {!$ok} break
		incr dbstart -2000
		incr dbend 2000
	}
	dict set dbobj prevdbloc $prevdbloc
	dict set dbobj fdbline $fdbline
	return $dblist
}

proc annotatemir_makegeneobj {genomef dbline} {
	global adata
	foreach {
		dchrom dstart dend geneobj strand mature1start mature1end loopstart loopend mature2start mature2end genename
	} $dbline break
	set temp [list_remove [lrange $dbline 5 10] {}]
	if {$temp ne [lsort -integer $temp]} {error "error in $dbline: values not in correct order"}
	if {$strand eq "-"} {set complement 1} else {set complement 0}
	set annotlist {}
	if {!$complement} {
		set temp {upstream 1}
	} else {
		set temp {downstream 1}
	}
	lappend temp [expr {$dstart-2000}] [expr {$dstart-100}]
	lappend annotlist $temp
	if {!$complement} {set side 5p} else {set side 3p}
	lappend annotlist [list flank$side 1 [expr {$dstart-100}] $dstart]
	if {$mature1start ne "" && $mature1end ne ""} {
		if {$dstart < $mature1start} {
			set code arm${side}
			if {$complement} {append code -mt}
			lappend annotlist [list $code 0 $dstart $mature1start]
		}
		if {$strand eq "+"} {
			lappend annotlist [list mature${side} 0 $mature1start [expr {$mature1start+1}]]
			lappend annotlist [list seed$side 0 [expr {$mature1start+1}] [expr {$mature1start+7}]]
			lappend annotlist [list mature${side}-sd 0[expr {$mature1start+7}] $mature1end]
		} else {
			lappend annotlist [list mature${side}-sd 0 $mature1start [expr {$mature1end-7}]]
			lappend annotlist [list seed$side 0 [expr {$mature1end-7}] [expr {$mature1end-1}]]
			lappend annotlist [list mature${side} 0 [expr {$mature1end-1}] $mature1end]
		}
		if {$mature1start < $mature1end} {
			set code arm${side}
			if {!$complement} {append code -mt}
			lappend annotlist [list $code 0 $mature1end $loopstart]
		}
	} else {
		lappend annotlist [list arm$side 0 $dstart $loopstart]
	}
	lappend annotlist [list loop 0 $loopstart $loopend]
	if {$complement} {set side 5p} else {set side 3p}
	if {$mature2start ne "" && $mature2end ne ""} {
		if {$dstart < $mature2start} {
			set code arm${side}
			if {$complement} {append code -mt}
			lappend annotlist [list $code 0 $loopend $mature2start]
		}
		if {$strand eq "+"} {
			lappend annotlist [list mature${side} 0 $mature2start [expr {$mature2start+1}]]
			lappend annotlist [list seed$side 0 [expr {$mature2start+1}] [expr {$mature2start+7}]]
			lappend annotlist [list mature${side}-sd 0[expr {$mature2start+7}] $mature2end]
		} else {
			lappend annotlist [list mature${side}-sd 0 $mature2start [expr {$mature2end-7}]]
			lappend annotlist [list seed$side 0 [expr {$mature2end-7}] [expr {$mature2end-1}]]
			lappend annotlist [list mature${side} 0 [expr {$mature2end-1}] $mature2end]
		}
		if {$mature2start < $mature2end} {
			set code arm${side}
			if {!$complement} {append code -mt}
			lappend annotlist [list $code 0 $mature2end $dend]
		}
	} else {
		lappend annotlist [list arm$side 0 $loopend $dend]
	}
	lappend annotlist [list flank$side 2 $dend [expr {$dend+100}]]
	if {$complement} {
		set temp {upstream 2}
	} else {
		set temp {downstream 2}
	}
	lappend temp [expr {$dend+100}] [expr {$dend+2000}]
	lappend annotlist $temp
	unset -nocomplain adata
	set adata(genomef) $genomef
	set adata(complement) $complement
	set adata(start) $dstart
	set adata(end) $dend
	set adata(genename) $genename
	set adata(chrom) $dchrom
	set adata(annotlist) $annotlist
puts -----${genename}-----
foreach a $adata(annotlist) {puts $a}
	return [array get adata]
}

proc open_mirfile {dbfile headerVar dpossVar {genecol name} {transcriptcol isomir}} {
	upvar $headerVar header
	upvar $dpossVar dposs
	set df [gzopen $dbfile]
	set header [tsv_open $df comment]
	set dposs [tsv_basicfields $header 3]
	set deffields {strand	mature1start mature1end	loopstart	loopend	mature2start mature2start}
	lappend deffields $genecol $transcriptcol
	lappend dposs {*}[list_cor $header $deffields]
	if {[lsearch [lrange $dposs 0 end-2] -1] != -1} {
		close $df
		puts stderr "error: gene file $dbfile misses the following fields: [list_sub $deffields [list_find [lrange $dposs 0 end-2] -1]]"
		exit 1
	}
	return $df
}

proc annotatemir {file genomefile dbfile name annotfile {genecol name} {transcriptcol isomir} {flanksize 100}} {
# putsvars file genomefile dbfile name annotfile genecol
	global genomef
	annot_init
	if {[catch {eof $genomef}]} {
		set genomef [genome_open $genomefile]
	}
	catch {close $f}; catch {close $df}; catch {close [dict get $dbobj df]} ;catch {close $o};
	set f [gzopen $file]
	set header [tsv_open $f comment]
	if {[catch {set poss [tsv_basicfields $header]}]} {
		if {[catch {set poss [tsv_basicfields $header 4]}]} {
			set poss [tsv_basicfields $header 3]
			lappend poss -1 -1 -1
		} else {
			lappend poss -1 -1
		}
		set noref 1
	} else {
		set noref 0
	}
	set fields [list_sub $header $poss]
	set dbobj [annotatemir_dbopen $dbfile $genecol $transcriptcol]
	set o [open $annotfile.temp w]
	puts -nonewline $o [join [list_fill [expr {[llength [split $comment \n]]-1}] \n]]
	set nh [list ${name}_impact ${name}_mir]
	puts $o [join $nh \t]
	set empty [join [list_fill [llength $nh] {}] \t]
	set dblist {}
	set counter 0
	set prevloc ""
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} {
			if {[eof $f]} break
			error "error: empty line in $file"
		}
		set loc [list_sub $line $poss]
		set ploc [lrange $loc 0 2]
		if {[lloc_compare $prevloc $ploc] > 0} {
			error "Cannot annotate because the variant file is not correctly sorted (sort correctly using \"cg select -s -\")"
		}
		set prevloc $ploc
		foreach {chr start end type ref alt} $loc break
		if {$type eq ""} {
			lset loc 3 del
			set type del
		}
		if {$start > $end} {
			puts stderr "location start > end error: $loc"
			exit 1
		}
		if {$noref} {
			switch $type {
				snp {
					set ref N; set alt N
					lset loc end N
				}
				default {set ref [expr {$end-$start}]}
			}
			lset loc end-1 $ref
		}
		incr counter
		if {$counter > 50000} {
			putslog $chr:$start
			set counter 0
		}
		# add all overlapping to dblist
		annotatemir_dbgetlist dbobj dblist $chr $start $end
		# join [list_subindex $dblist 4] \n\n
		# check for multiple alleles, process these separately (alist contains >1 loc)
		set malt [split $alt ,]
		if {[llength $malt] > 1} {
			set alist {}
			foreach a $malt {
				lappend alist [lreplace $loc end end $a]
			}
		} else {
			set alist [list $loc]
		}
		set ahitgenes {}
		# separately for each alt allele
		# check for overlap, mark genes from dblist for removal 
		# that are before current var, annotate
		foreach loc $alist {
			set num 0
			set remove {}
			set hitgenes ""
			# join [list_subindex $dblist 4] \n\n
			foreach dbline $dblist {
				foreach {dc ds de geneobj} $dbline break
				incr ds -2000
				incr de 2000
				set chrcompar [chr_compare $dc $chr]
				if {$chrcompar < 0} {
					lappend remove $num
				} elseif {$chrcompar == 0} {
					if {$de <= $start} {
						lappend remove $num
					} elseif {$ds < $end} {
						if {[catch {dict get $geneobj end}]} {
							set geneobj [annotatemir_makegeneobj $genomef $dbline]
							lset dblist $num 3 $geneobj
						}
						set result [annotatemir_one $loc $geneobj]
						lappend hitgenes $result
					}
				}
				incr num
			}
			set result {}
			set impacts [list_subindex $hitgenes 0]
			if {[llength $hitgenes] == 1} {
				set result [lindex $hitgenes 0]
			} elseif {[llength $hitgenes]} {
				set pos 0
				foreach h $nh {
					set list [list_subindex $hitgenes $pos]
					set nrlist [list_remdup $list]
					if {[llength $nrlist] == 1} {
						lappend result [lindex $nrlist 0]
					} else {
						lappend result [join $list \;]
					}
					incr pos
				}
			}
			lappend ahitgenes $result
		}
		# remove genes that have been passed by
		if {[llength $remove]} {
			set dblist [list_sub $dblist -exclude $remove]
		}
		# join results to one line
		if {[llength $ahitgenes] == 1} {
			set line [lindex $ahitgenes 0]
			if {[llength $line]} {
				set result [join $line \t]
			} else {
				set result $empty
			}
		} elseif {[llength $ahitgenes]} {
			set result {}
			set pos 0
			foreach h $nh {
				set list [list_subindex $ahitgenes $pos]
				set nrlist [list_remdup $list]
				if {[llength $nrlist] == 1} {
					lappend result [lindex $nrlist 0]
				} else {
					lappend result [join $list \,]
				}
				incr pos
			}
			set result [join $result \t]
		} else {
			set result $empty
		}
		puts $o $result
	}
	close $genomef
	close $o; catch {close $f};	catch {close $df}
	file rename -force $annotfile.temp $annotfile
}
