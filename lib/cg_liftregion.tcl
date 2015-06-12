proc cg_liftregion {args} {
#	set pos 0
#	set sumfields {}
#	foreach {key value} $args {
#		switch -- $key {
#			-sumfields {
#				set sumfields $value
#			}
#			-- break
#			default {
#				break
#			}
#		}
#		incr pos 2
#	}
#	set args [lrange $args $pos end]
	if {([llength $args] < 3)} {
		errorformat liftregion
		exit 1
	}
	foreach {file resultfile liftoverfile} $args break
	set unmappedfile $resultfile.unmapped
	set liftoverfile [liftoverfile $liftoverfile]
	if {[file exists $resultfile]} {
		error "file $resultfile already exists"
	}
	catch {gzclose $f} ; catch {gzclose $fl} ; catch {close $o} ; catch {close $ou}
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set chrpos [lindex $poss 0]
	set strandpos [lsearch $header strand]
	if {$strandpos != -1} {
		lappend poss $strandpos
	}
	set strand {}
	set fl [gzopen $liftoverfile]
	set lheader [tsv_open $fl]
	set lposs [list_cor $lheader {chromosome begin end strand destchromosome destbegin destend deststrand}]
	if {"-1" in $lposs} {
		error "liftoverfile $liftoverfile is missing following column(s): [list_sub {chromosome begin end strand destchromosome destbegin destend deststrand} [list_find $lposs -1]]"
	}
	set newheader [list_union [list_sub $header $poss] $header]
	lappend newheader beforeliftover
	set tempfile [tempfile]
	set o [open $tempfile w]
	puts $o "# liftregion from $file"
	puts $o "# using $liftoverfile"
	puts $o [join $newheader \t]
	set ou [open $unmappedfile.temp w]
	puts $ou "# unmapped liftregion from $file"
	puts $ou "# using $liftoverfile"
	puts $ou [join [list {*}[list_union [list_sub $header $poss] $header] oriregion] \t]
	set liftregions {}
	set lline [list_sub [split [gets $fl] \t] $lposs]
	lset lline 4 [chr_clip [lindex $lline 4]]
	foreach {srcchromosome srcbegin srcend srcstrand destchromosome destbegin destend deststrand} $lline break
	while 1 {
		if {[gets $f oline] == -1} break
		set line [split $oline \t]
		set loc [list_sub $line $poss]
		foreach {chromosome begin end strand} $loc break
		lset line $chrpos [chr_clip $chromosome]
		set restline [list_sub $line -exclude $poss]
		set before ${chromosome}-${begin}-${end}
		if {$strand ne ""} {append before -$strand}
		lappend restline $before
		set restline [join $restline \t]
		set pos 0
		set comp [reg_compare $lline $loc]
		if {$comp <= 0} {
			while {$lline ne ""} {
				if {$comp == 0} {lappend liftregions $lline}
				if {[gets $fl lline] == -1} break
				set lline [list_sub [split $lline \t] $lposs]
				lset lline 4 [chr_clip [lindex $lline 4]]
				set comp [reg_compare $lline $loc]
				if {$comp > 0} {
					break
				}
			}
		}
		set pos 0
		foreach templline $liftregions {
			set comp [reg_compare $templline [list $chromosome $begin $end]]
			# putsvars templline begin end
			if {$comp < 0} {
				incr pos
				continue
			} elseif {$comp == 0} {
				# overlap
				foreach {srcchromosome srcbegin srcend srcstrand destchromosome destbegin destend deststrand} $templline break
				if {$begin < $srcbegin} {
					if {$strandpos != -1} {
						puts $ou $chromosome\t$begin\t$srcbegin\t$strand\t$restline
					} else {
						puts $ou $chromosome\t$begin\t$srcbegin\t$restline
					}
					set begin $srcbegin
				}
				if {$end <= $srcend} {
					set ubegin $begin ; set uend $end
					set begin $end
				} elseif {$end > $srcend} {
					set ubegin $begin ; set uend $srcend
					set begin $srcend
				}
				# putsvars ubegin uend
				if {$uend > $ubegin} {
					if {$srcstrand eq $deststrand} {
						set ustrand $strand
						set ubegin [expr {$ubegin + $destbegin - $srcbegin}]
						set uend [expr {$uend + $destbegin - $srcbegin}]
					} else {
						if {$strand eq "+"} {set ustrand "-"} else {set ustrand "+"}
						set keepuend $uend
						set uend [expr {$destend - $ubegin + $srcbegin}]
						set ubegin [expr {$destend - $keepuend + $srcbegin}]
					}
					if {$strandpos != -1} {
						puts $o $destchromosome\t$ubegin\t$uend\t$ustrand\t$restline
					} else {
						puts $o $destchromosome\t$ubegin\t$uend\t$restline
					}
				}
			}
		}
		if {$end > $begin} {
			if {$strandpos != -1} {
				puts $ou $chromosome\t$begin\t$end\t$strand\t$restline
			} else {
				puts $ou $chromosome\t$begin\t$end\t$restline
			}
		}
		set liftregions [lrange $liftregions $pos end]
	}
	gzclose $f ; gzclose $fl ; close $o ; close $ou
	#
	# sort result
	#
	set tempfile2 [tempfile]
	cg select -s {chromosome begin end beforeliftover} $tempfile $tempfile2
	# collapse overlapping regions
	set o [open $resultfile.temp w]
	puts $o "# liftregion from $file"
	puts $o "# using $liftoverfile"
	close $o
	cg regcollapse $tempfile2 >> $resultfile.temp
	#
	# rename result, cleanup
	#
	file rename -force $resultfile.temp $resultfile
	file delete $tempfile $tempfile2
	file rename -force $unmappedfile.temp $unmappedfile
}
