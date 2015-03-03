#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc annotatemir_one {oneloc geneobj {addtranscriptname 1}} {
# putsvars oneloc geneobj
	foreach {chrom snpstart snpend} $oneloc break
	if {$snpstart == $snpend} {
		incr snpstart -1
		incr snpend
	}
	unset -nocomplain adata
	array set adata $geneobj
	set complement $adata(complement)
	set dbstart [lindex $adata(annotlist) 0 3]
	set dbend [lindex $adata(annotlist) end 4]
	set hstart $adata(start)
	set hend $adata(end)
	set flank1 $adata(flank1)
	set flank2 $adata(flank2)
	if {$flank2 eq ""} {set flank2 $flank1}
	array set da {}
	if {$snpstart <= $hstart && $snpend >= $hend} {
		set impact GENEDEL
	} else {
		set impact {}
		set dist {}
		set snpsize [expr {$snpend - $snpstart - 1}]
		list_foreach {annot ref distback b e} $adata(annotlist) {
			if {$snpend > $b && $snpstart < $e} {
				if {$snpend >= $e} {
					set tsnpend $e
					set endend 1
				} else {
					set tsnpend $snpend
					set endend 0
				}
				if {$snpstart <= $b} {
					set tsnpstart $b
					set startend 1
				} else {
					set tsnpstart $snpstart
					set startend 0
				}
				set temp $annot
				set snpto {}
				if {$distback == -1} {
					if {!$complement} {set sign -} else {set sign +}
					set num1 [expr {$e - $tsnpstart}]
					set num2 [expr {$e - $tsnpend + 1}]
				} elseif {$distback == +1} {
					if {!$complement} {set sign +} else {set sign -}
					set num1 [expr {$tsnpstart - $b + 1}]
					set num2 [expr {$tsnpend - $b}]
				} elseif {!$complement} {
					set sign +
					set num1 [expr {$tsnpstart - $b + 1}]
					set num2 [expr {$tsnpend - $b}]
				} else {
					set sign +
					set num1 [expr {$e - $tsnpstart}]
					set num2 [expr {$e - $tsnpend + 1}]
				}
				if {$annot eq "upstream"} {
					if {!$complement} {
						incr num1 $flank1
						incr num2 $flank1
					} else {
						incr num1 $flank2
						incr num2 $flank2
					}
				} elseif {$annot eq "downstream"} {
					if {!$complement} {
						incr num1 $flank2
						incr num2 $flank2
					} else {
						incr num1 $flank1
						incr num2 $flank1
					}
				}
				set unum1 $num1
				set unum2 $num2
				if {$num1 == $num2} {
# putsvars snpstart snpend annot ref distback b e num1 num2 startend endend complement
					if {$startend} {
						if {!$complement} {set unum1 e$unum1} else {set unum1 ${unum1}e}
					}
					if {$endend} {
						if {!$complement} {set unum1 ${unum1}e} else {set unum1 e${unum1}}
					}
					set temp $annot\($ref$sign$unum1\)
				} else {
					if {$startend} {
						if {!$complement} {set unum1 e$unum1} else {set unum1 ${unum1}e}
					}
					if {$endend} {
						if {!$complement} {set unum2 ${unum2}e} else {set unum2 e${unum2}}
					}
					if {!$complement} {
						set temp $annot\($ref$sign$unum1:$unum2\)
					} else {
						set temp $annot\($ref$sign$unum2:$unum1\)
					}
				}
# putsvars snpstart snpend b e complement num1 num2 temp
				if {[regexp ^mature $annot] && $num2 >= 2 && $num1 <= 7} {
					append temp seed
				}
				lappend impact $temp
			}
		}
	}
	if {$complement} {set impact [list_reverse $impact]}
	if {$addtranscriptname} {
		return [list $adata(transcriptname):[join $impact &] $adata(genename)]
	} else {
		return [list [join $impact &] $adata(genename)]
	}
}

proc annotatemir_dbopen {dbfile genecol transcriptcol} {
	set df [gzopen $dbfile]
	set header [tsv_open $df comment]
	set dposs [tsv_basicfields $header 3]
	set deffields {strand	mature1start mature1end	loopstart	loopend	mature2start mature2end}
	lappend deffields $genecol $transcriptcol
	lappend dposs {*}[list_cor $header $deffields]
	if {[lsearch [lrange $dposs 0 end-2] -1] != -1} {
		close $df
		error "error: gene file $dbfile misses the following fields: [list_sub $deffields [list_find [lrange $dposs 3 end-2] -1]]"
	}
	# insert empty space for genobj
	lappend dposs -1
	set geneobjpos [llength $dposs]
	incr geneobjpos -1
	while {![eof $df]} {
		set fdbline [split [gets $df] \t]
		set fdbline [list_sub $fdbline $dposs]
		foreach {dbchr dbstart dbend} $fdbline break
		if {![isint $dbstart] || ![isint $dbend]} continue
		break
	}
	set dbobj [dict create df $df dposs $dposs fdbline $fdbline prevdbloc [lrange $fdbline 0 2] geneobjpos $geneobjpos]
	return $dbobj
}

proc annotatemir_dbgetlist {dbobjVar dblistVar chr start end} {
	upvar $dbobjVar dbobj
	upvar $dblistVar dblist
	set df [dict get $dbobj df]
	if {[eof $df]} {return $dblist}
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

proc annotatemir_makegeneobj {genomef dbline {flanksizes 100} {isomirname 0}} {
# putsvars genomef dbline flanksizes
	global adata
	set flank1 [lindex $flanksizes 0]
	set flank2 [lindex $flanksizes 1]
	if {$flank2 eq ""} {set flank2 $flank1}
	foreach {
		dchrom dstart dend strand mature1start mature1end loopstart loopend mature2start mature2end genename transcriptname
	} $dbline break
	set temp [list_remove [lrange $dbline 4 9] {}]
	if {$temp ne [lsort -integer $temp]} {error "error in $dbline: values not in correct order"}
	if {$strand eq "-"} {set complement 1} else {set complement 0}
	set annotlist {}
	if {!$complement} {
		set temp {upstream a -1}
	} else {
		set temp {downstream a -1}
	}
	lappend temp [expr {$dstart-2000}] [expr {$dstart-$flank1}]
	lappend annotlist $temp
	if {!$complement} {set side 5p} else {set side 3p}
	lappend annotlist [list flank$side a -1 [expr {$dstart-$flank1}] $dstart]
	if {$mature1start ne "" && $mature1end ne ""} {
		if {$dstart < $mature1start} {
			lappend annotlist [list arm${side} m -1 $dstart $mature1start]
		}
		if {$isomirname} {
			if {!$complement} {
				set isomir [expr {$mature1start-$dstart+1}]_[expr {$mature1end-$dstart+1}]
			} else {
				set isomir [expr {$dend-$mature1end+1}]_[expr {$dend-$mature1start+1}]
			}
		} else {
			set isomir {}
		}
		lappend annotlist [list mature${side}${isomir} a 0 $mature1start $mature1end]
		if {$mature1end < $loopstart} {
			lappend annotlist [list arm${side} m +1 $mature1end $loopstart]
		}
	} else {
		lappend annotlist [list arm$side l -1 $dstart $loopstart]
	}
	lappend annotlist [list loop a 0 $loopstart $loopend]
	if {$complement} {set side 5p} else {set side 3p}
	if {$mature2start ne "" && $mature2end ne ""} {
		if {$loopend < $mature2start} {
			lappend annotlist [list arm${side} m -1 $loopend $mature2start]
		}
		if {$isomirname} {
			if {!$complement} {
				set isomir [expr {$mature2start-$dstart+1}]_[expr {$mature2end-$dstart+1}]
			} else {
				set isomir [expr {$dend-$mature2end+1}]_[expr {$dend-$mature2start+1}]
			}
		} else {
			set isomir {}
		}
		lappend annotlist [list mature${side}${isomir} a 0 $mature2start $mature2end]
		if {$mature2end < $dend} {
			lappend annotlist [list arm${side} m +1 $mature2end $dend]
		}
	} else {
		lappend annotlist [list arm$side l +1 $loopend $dend]
	}
	lappend annotlist [list flank$side a +1 $dend [expr {$dend+$flank2}]]
	if {$complement} {
		set temp {upstream a +1}
	} else {
		set temp {downstream a +1}
	}
	lappend temp [expr {$dend+$flank2}] [expr {$dend+2000}]
	lappend annotlist $temp
	unset -nocomplain adata
	set adata(genomef) $genomef
	set adata(complement) $complement
	set adata(start) $dstart
	set adata(end) $dend
	set adata(genename) $genename
	if {$transcriptname ne ""} {
		set adata(transcriptname) $transcriptname
	} else {
		set adata(transcriptname) $genename
	}
	set adata(chrom) $dchrom
	set adata(annotlist) $annotlist
	set adata(flank1) $flank1
	set adata(flank2) $flank2
	# puts -----${genename}-----
	# foreach a $adata(annotlist) {puts $a}
	return [array get adata]
}

# about the mirvas parameter:
# for genomecomb the parameter must must be 0: 
# for mirvas it will be 1 (table only) or 2 (include graphical results as well)
# Although this piece of code is shared with mirvas, mirvas does several things differently.
# genomecomb does not have all the tools to do the structural analysis.
# 
proc annotatemir {file genomefile dbfile name resultfile {genecol name} {transcriptcol transcript} {flanksizes 100} {isomirname 0} {mirvas 0}} {
# putsvars file genomefile dbfile name resultfile genecol transcriptcol flanksizes mirvas
	global genomef
	set file [file normalize $file]
	set genomefile [file normalize $genomefile]
	set dbfile [file normalize $dbfile]
	set resultfile [file normalize $resultfile]
	file mkdir [file dir $resultfile]
	if {!$mirvas} {annot_init}
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
	set geneobjpos [dict get $dbobj geneobjpos]
	set o [open $resultfile.temp w]
	puts -nonewline $o $comment
	if {$mirvas} {
		set addtranscriptname 0
		set nh {chromosome begin end type ref alt mir_location mir_name}
		set temp [annotatemir_one_struct_fields]
		lappend nh {*}$temp
		lappend nh {*}[list_sub $header -exclude $poss]
		set mirvasempty [list_fill [llength $temp] {}]
		set mo [mirvas_start $mirvas $file $resultfile $flanksizes]
	} else {
		set addtranscriptname 1
		set nh [list ${name}_impact ${name}_mir]
	}
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
		if {$mirvas} {
			progress_step 1.0
			set restline [list_sub $line -exclude $poss]
		}
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
		if {$mirvas && [llength $alist] > 1} {progress_step [expr {1.0/[llength $alist]}]}
		foreach loc $alist {
			log "annotating $loc"
			set num 0
			set remove {}
			set hitgenes ""
			# join [list_subindex $dblist 4] \n\n
			set found 0
			set tododblist {}
			foreach dbline $dblist {
				foreach {dc ds de} $dbline break
				incr ds -2000
				incr de 2000
				set chrcompar [chr_compare $dc $chr]
				if {$chrcompar < 0} {
					lappend remove $num
				} elseif {$chrcompar == 0} {
					if {$de <= $start} {
						lappend remove $num
					} elseif {$ds < $end} {
						set geneobj [lindex $dbline $geneobjpos]
						if {[catch {dict get $geneobj end}]} {
							set geneobj [annotatemir_makegeneobj $genomef $dbline $flanksizes $isomirname]
							lset dbline $geneobjpos $geneobj
							lset dblist $num $dbline
						}
						lappend tododblist $dbline
						set found 1
					}
				}
				incr num
			}
			if {$mirvas && [llength $tododblist] > 1} {
				progress_step [expr {[progress_step]/[llength $tododblist]}]
			}
			foreach dbline $tododblist {
				set geneobj [lindex $dbline $geneobjpos]
				set result [annotatemir_one $loc $geneobj $addtranscriptname]
				if {$mirvas} {
					# putsvars loc geneobj resultfile mo mirvas result
					lappend result {*}[annotatemir_one_struct $loc $geneobj $resultfile $mo $mirvas $result]
					puts $o [join [list_concat $loc $result $restline] \t]
				}
				lappend hitgenes $result
			}
			if {$mirvas} {
				if {!$found} {puts $o [join [list_concat $loc {{} {}} $mirvasempty $restline] \t]}
				continue
			}
			set result {}
			set impacts [list_subindex $hitgenes 0]
			set hitgenes [list_remdup $hitgenes]
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
		if {$mirvas} continue
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
	# genome_close $genomef
	if {$mirvas == 2} {
		mirvas_draw_close $mo
	}
	close $o; catch {close $f};	catch {close $df}
	file rename -force $resultfile.temp $resultfile
}
