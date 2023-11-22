proc multigene_findfield {header tryfields msg fieldnameVar {file {}}} {
	upvar $fieldnameVar fieldname
	foreach fieldname $tryfields {
		set pos [lsearch $header $fieldname]
		if {$pos != -1} break
	}
	if {$pos == -1} {
		set error "Could not find $msg field (checked for: $tryfields)"
		if {$file ne ""} {append error " in file $file"}
		error $error
	}
	return $pos
}

proc multigene_genecounts_moutput {merged files} {
	upvar gcounta gcounta
	upvar dummya dummya
	set genes [list_subindex $merged 4]
	set line [lindex $merged 0]
	lset line 2 [lindex $merged end 2]
	set num -1
	foreach file $files {
		incr num
		set found 0
		foreach gene $genes {
			if {[info exists gcounta($num,$gene)]} {
				lappend line {*}$gcounta($num,$gene)
				set found 1
				break
			}
		}
		if (!$found) {
			lappend line {*}$dummya($num)
		}
	}
	return $line
}

proc cg_multigene {args} {
	set novelpattern novel
	cg_options multigene args {
		-novelpattern - -match {
			set novelpattern $value
		}
	} genecounts 2
	analysisinfo_combine $genecounts $args
	unset -nocomplain geneinfoa
	unset -nocomplain gcounta
	unset -nocomplain dummya
	unset -nocomplain chrgenea
	unset -nocomplain novela
	set f [gzopen [lindex $args 0]]
	set header [tsv_open $f comments1]
	gzclose $f
	set poss [list_sub [tsv_basicfields $header 7 0] {0 1 2 6}]
	set common {chromosome begin end strand}
	set basicfieldspresent [list_sub $common -exclude [list_find $poss -1]]
	set newheader {gene geneid}
	set num -1
	foreach file $args {
		incr num
		set f [gzopen $file]
		set header [tsv_open $f comment]
		set poss [list_sub [tsv_basicfields $header 7 0] {0 1 2 6}]
		set chrpos [lindex $poss 0]
		if {[llength $basicfieldspresent] < 4} {
			list_addnew basicfieldspresent {*}[list_sub $common -exclude [list_find $poss -1]]
		}
		set genepos [multigene_findfield $header {gene_name genename gene name2 geneid gene_id} gene fieldname $file]
		set geneidpos [multigene_findfield $header {gene_id geneid gene} geneid fieldname $file]
		lappend poss $genepos $geneidpos
		set restposs [list_find -regexp $header -]
		if {[llength $restposs] > 0} {
			set restfields [list_sub $header $restposs]
			lappend newheader {*}$restfields
		} else {
			set restfields [list_sub $header -exclude [list_remove [list_remdup $poss] -1]]
			set restposs [list_cor $header $restfields]
			foreach field $restfields {
				lappend newheader ${field}-[file_rootname $file]
			}
		}
		set dummya($num) {}
		foreach field $restfields {lappend dummya($num) 0}
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			foreach {chr begin end strand} [list_sub $line $poss] break
			if {$chrpos != -1} {
				set chr [chr_clip $chr]
				lset line $chrpos $chr
			}
			set gene [lindex $line $genepos]
			if {$novelpattern ne "" && [regexp $novelpattern $gene]} {
				set gene [gene_name $chr $strand $begin $end]
				if {![info exists geneinfoa($gene)]} {
					if {[info exists novela($chr,$strand)]} {
						foreach {b e g} $novela($chr,$strand) {
							if {$begin < $e && $end >= $b} {
								set gene $g
								break
							}
						}
						if {![info exists geneinfoa($gene)]} {
							lappend novela($chr,$strand) $begin $end $gene
						}
					} else {
						lappend novela($chr,$strand) $begin $end $gene
					}
				} else {
					lappend novela($chr,$strand) $begin $end $gene
				}
			}
			if {[info exists chrgene($gene-$chr)]} {
				set gene gene-$chr
			} elseif {[info exists chrgenea($gene)]} {
				if {$chrgenea($gene) ne $chr} {
					set gene gene-$chr
					set chrgenea($gene) $chr
				}
			} else {
				set chrgenea($gene) $chr
			}
			list_addnew geneinfoa($gene) [list_sub $line $poss]
			set gcounta($num,$gene) [list_sub $line $restposs]
		}
		gzclose $f
	}
	# parray geneinfoa
	# parray gcounta
	catch {close $o}
	set o [wgzopen $genecounts.temp[gzext $genecounts]]
	set basicfieldspresent [list_sub $basicfieldspresent [list_remove [list_cor $basicfieldspresent $common] -1]]
	set poss [list_cor $common $basicfieldspresent]
	if {[llength $basicfieldspresent]} {
		set temp [join $basicfieldspresent \t]\t
	} else {
		set temp {}
	}
	puts $o $comments1$temp[join $newheader \t]
	unset -nocomplain merge
	set genes [bsort [array names geneinfoa]]
	foreach gene $genes {
		set lines $geneinfoa($gene)
		set chr [list_remdup [list_subindex $lines 0]]
		if {$chr eq {{}}} {
			# no location
		} elseif {[llength $lines] > 1} {
			if {[llength $chr] > 1} {
				# should never get here -> solved earlier using chrgenea
				puts stderr "error combining gene count files: gene $gene is in different chromosomes ($chr)"
			}
			set begin [list_remove [list_subindex $lines 1] {}]
			if {[llength $begin] > 1} {
				set begin [lmath_min $begin]
			}
			set end [list_remove [list_subindex $lines 2] {}]
			if {[llength $end] > 1} {
				set end [lmath_max $end]
			}
			set temp [list_remove [list_subindex $lines 3] {}]
			set strand [join [list_remdup [split $temp { ,}]] ,]
			set line [lindex $lines 0]
			lset line 1 $begin
			lset line 2 $end
			lset line 3 $strand
			if {$novelpattern ne "" && [regexp $novelpattern $gene]} {
				set ogene [gene_name $chr $strand $begin $end]
				lset line 4 $ogene
				lset line 5 $ogene
			}
			set lines [list $line]
		} else {
			set line [lindex $lines 0]
			if {$novelpattern ne "" && [regexp $novelpattern $gene]} {
				foreach {chr begin end strand} $line break
				set ogene [gene_name $chr $strand $begin $end]
				lset line 4 $ogene
				lset line 5 $ogene
			}
			set lines [list $line]
		}
		foreach line $lines {
			set geneid [lindex $line end]
			if {$gene in ". dummy"} continue
			set num -1
			foreach file $args {
				incr num
				if {[info exists gcounta($num,$gene)]} {
					lappend line {*}$gcounta($num,$gene)
				} else {
					lappend line {*}$dummya($num)
				}
			}
			if {[llength $basicfieldspresent] < 4} {
				set line [list {*}[list_sub $line $poss] {*}[lrange $line 4 end]]
			}
			puts $o [join $line \t]
		}
	}
	close $o
	cg select -overwrite 1 -s - $genecounts.temp[gzext $genecounts] $genecounts.temp2[gzext $genecounts]
	file delete $genecounts.temp[gzext $genecounts]
	file rename -force $genecounts.temp2[gzext $genecounts] $genecounts
}
