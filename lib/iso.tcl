proc exons_startsends2list {starts ends {sizeVar {}}} {
	if {$sizeVar ne ""} {upvar $sizeVar size}
	set size 0
	set exons {}
	foreach begin [split $starts ,] end [split $ends ,] {
		if {$begin eq ""} continue
		set size [expr {$size + ($end - $begin)}]
		lappend exons $begin-$end
	}
	return $exons
}

proc iso_name {chromosome strand exonStarts exonEnds {sizeVar {}}} {
	if {$sizeVar ne ""} {upvar $sizeVar size}
	set size 0
	set starts [split [string trimright $exonStarts ,] ,]
	set begin [lindex $starts 0]
	set newname novelt_${chromosome}_${begin}${strand}
	set pe -1
	foreach s $starts e [split [string trimright $exonEnds ,] ,] {
		if {$s eq ""} continue
		if {$pe != -1} {
			append newname i[expr {$s-$pe}]
		}
		set pe $e
		set esize [expr {$e-$s}]
		append newname e$esize
		set size [expr {$size + $esize}]
	}
	return $newname
}

proc iso_combine_findfield {header tryfields msg fieldnameVar {file {}}} {
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

proc iso_combine_genecounts_moutput {merged files} {
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

proc cg_iso_combine_genecounts {genecounts args} {
	analysisinfo_combine $genecounts $args
	unset -nocomplain geneinfoa
	unset -nocomplain gcounta
	unset -nocomplain dummya
	unset -nocomplain chrgenea
	set f [gzopen [lindex $args 0]]
	set header [tsv_open $f comments1]
	close $f
	set poss [list_sub [tsv_basicfields $header 7 0] {0 1 2 6}]
	set common [list_sub $header $poss]
	set newheader $common
	lappend newheader gene geneid
	set num -1
	foreach file $args {
		incr num
		set f [gzopen $file]
		set header [tsv_open $f comment]
		set poss [list_cor $header $common]
		set genepos [iso_combine_findfield $header {gene_name gene name2 geneid gene_id} gene fieldname $file]
		set geneidpos [iso_combine_findfield $header {gene_id geneid gene} geneid fieldname $file]
		lappend poss $genepos $geneidpos
		set restposs [list_find -regexp $header -]
		set restfields [list_sub $header $restposs]
		lappend newheader {*}$restfields
		set dummya($num) {}
		foreach field $restfields {lappend dummya($num) 0}
		while {[gets $f line] != -1} {
			set line [split $line \t]
			set chr [lindex $line 0]
			set gene [lindex $line $genepos]
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
		close $f
	}
	set o [open $genecounts.temp w]
	puts $o $comments1[join $newheader \t]
	set merge(+) {}
	set merge(-) {}
	set merge(.) {}
	set merge() {}
	set genes [bsort [array names geneinfoa]]
	foreach gene $genes {
		set lines $geneinfoa($gene)
		if {[llength $lines] > 1} {
			set chr [list_remdup [list_subindex $lines 0]]
			if {[llength $chr] > 1} {
				# should never get here -> solved earlier using chrgenea
				puts stderr "error combining gene count files: gene $gene is in different chromosomes ($chr)"
			}
			set begin [lmath_min [list_subindex $lines 1]]
			set end [lmath_max [list_subindex $lines 2]]
			set strand [join [list_remdup [split [list_subindex $lines 3] { ,}]] ,]
			set line [lindex $lines 0]
			lset line 1 $begin
			lset line 2 $end
			lset line 3 $strand
		} else {
			set line [lindex $lines 0]
		}
		set geneid [lindex $line end]
		if {$gene in ". dummy"} continue
		if {[regexp novel $gene]} {
			foreach {chr begin end strand} $line break
			if {![llength $merge($strand)]} {
				lappend merge($strand) $line
				continue
			} else {
				set ready 0
				foreach {prevchr prevbegin} [lindex $merge($strand) 0] break
				set prevend [lindex $merge($strand) end 2]
				if {$chr ne $prevchr || $begin >= $prevend} {
					set merged $merge($strand)
					set merge($strand) [list $line]
					set line [iso_combine_genecounts_moutput $merged $args]
				} else {
					lappend merge($strand) $line
					continue
				}
			}
		} else {
			set num -1
			foreach file $args {
				incr num
				if {[info exists gcounta($num,$gene)]} {
					lappend line {*}$gcounta($num,$gene)
				} else {
					lappend line {*}$dummya($num)
				}
			}
		}
		puts $o [join $line \t]
	}
	if {[llength $merge(+)]} {
		puts $o [join [iso_combine_genecounts_moutput $merge(+) $args] \t]
	}
	if {[llength $merge(-)]} {
		puts $o [join [iso_combine_genecounts_moutput $merge(-) $args] \t]
	}
	close $o
	cg select -overwrite 1 -s - $genecounts.temp $genecounts.temp2
	file delete $genecounts.temp
	file rename -force $genecounts.temp2 $genecounts
}

proc iso_combine_job {projectdir isocaller {iso_match {}}} {
	upvar job_logdir job_logdir
	# combined analysis
	cd $projectdir
	mkdir compar
	set exproot [file tail $projectdir]
	if {$isocaller eq "*"} {
		set root $exproot
	} else {
		set root ${isocaller}-$exproot
	}
	set isoformfiles [bsort [jobglob samples/*/isoform_counts-${isocaller}-*.tsv]]
	if {[llength $isoformfiles]} {
		job iso_compar-isoform_counts-$root \
		-deps $isoformfiles \
		-targets {
			compar/isoform_counts-$root.tsv
		} -vars {
			isoformfiles exproot root isocaller iso_match
		} -code {
			set isoformcounts compar/isoform_counts-$root.tsv
			cg multitranscript -match $iso_match $isoformcounts {*}$isoformfiles
		}
	}
	set genefiles [bsort [jobglob samples/*/gene_counts-${isocaller}-*.tsv]]
	if {[llength $genefiles]} {
		job iso_compar-gene_counts-$root \
		-deps $genefiles \
		-targets {
			compar/gene_counts-$root.tsv
		} -vars {
			genefiles exproot root isocaller
		} -code {
			set genecounts compar/gene_counts-$root.tsv
			cg_iso_combine_genecounts $genecounts {*}$genefiles
		}
	}
	set totalcountsfiles [bsort [jobglob samples/*/totalcounts-${isocaller}-*.tsv]]
	if {[llength $totalcountsfiles]} {
		job iso_compar-totalcounts-$root \
		-deps $totalcountsfiles \
		-targets {
			compar/totalcounts-$root.tsv
		} -vars {
			totalcountsfiles exproot root isocaller
		} -code {
			analysisinfo_combine compar/totalcounts-$root.tsv $totalcountsfiles
			cg paste {*}$totalcountsfiles > compar/totalcounts-$root.tsv
		}
	}
}
