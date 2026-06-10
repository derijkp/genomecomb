proc checkexonoverlaps {overlaps exonStarts exonEnds} {
	set maxeoverlap 0
	set newoverlaps {}
	set exonStarts [split [string trim $exonStarts ,] ,]
	set exonEnds [split [string trim $exonEnds ,] ,]
	foreach line $overlaps {
		foreach {tgeneid tgene overlap ipct opct obegin oend texonStarts texonEnds} $line break
		set texonStarts [split [string trim $texonStarts ,] ,]
		set texonEnds [split [string trim $texonEnds ,] ,]
		# not most efficient, but will do for now
		set eoverlap 0
		# putsvars exonStarts exonEnds texonStarts texonEnds
		foreach b $exonStarts e $exonEnds {
			foreach tb $texonStarts te $texonEnds {
				if {$tb >= $e} break
				if {$te < $b} continue
				set ob [max $b $tb]
				set oe [min $e $te]
				set eoverlap [expr {$eoverlap + $oe - $ob}]
			}
		}
		lappend line $eoverlap
		if {$eoverlap == $maxeoverlap} {
			lappend newoverlaps $line
		} elseif {$eoverlap > $maxeoverlap} {
			set maxeoverlap $eoverlap
			set newoverlaps [list $line]
		}
	}
	if {$maxeoverlap > 0} {
		return $newoverlaps
	} else {
		return $overlaps
	}
}

proc cg_fix_isogene {args} {
	set file {}
	set genefile {}
	set renamefile {}
	cg_options fix_isogene args {
		-genefile {set genefile $value}
		-renamefile {set renamefile $value}
	}  {file outfile} 0 2
	if {[llength $args] > 2} {
		errorformat fix_isogene
	}
	if {$genefile eq "" && $renamefile eq ""} {
		error "error: either option -genefile or -renamefile must be given"
	}
	if {![file exists $file]} {
		error "error fixing \"$file\": file does not exist"
	}

	if {$genefile ne ""} {
		# load gene data
		catch {close $f}
		set f [gzopen $genefile]
		set header [tsv_open $f]
		set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 12 13 7 8}]
		# list_sub $header $poss
		unset -nocomplain chra
		unset -nocomplain begina
		unset -nocomplain enda
		unset -nocomplain stranda
		unset -nocomplain genea
		unset -nocomplain exonStartsa
		unset -nocomplain exonEndsa
		unset -nocomplain chrstra
		while {[gets $f line] != -1} {
			foreach {chr begin end strand gene geneid exonStarts exonEnds} [list_sub [split $line \t] $poss] break
			# putsvars chr begin end strand gene geneid
			if {![info exists chra($geneid)]} {
				set begina($geneid) $begin
				set enda($geneid) $end
				set stranda($geneid) $strand
				set genea($geneid) $gene
				set exonStartsa($geneid) $exonStarts
				set exonEndsa($geneid) $exonEnds
				set chra($geneid) $chr
				lappend chrstra($chr,$stranda($geneid)) $geneid
			} else {
				if {$begin < $begina($geneid)} {set begina($geneid) $begin}
				if {$end > $enda($geneid)} {set enda($geneid) $end}
				if {$strand != $stranda($geneid)} {error "different strand for $geneid: $line"}
				if {$chr != $chra($geneid) && [chr_clip $chr] ne "Y"} {error "different chr for $geneid: $line"}
				if {$gene != $genea($geneid)} {error "different gene for $geneid: $line"}
			}
		}
		gzclose $f
		unset -nocomplain checka
		foreach chrstr [array names chrstra] {
			set list {}
			foreach geneid $chrstra($chrstr) {
				lappend list [list $begina($geneid) $enda($geneid) $geneid $genea($geneid) $exonStartsa($geneid) $exonEndsa($geneid)]
			}
			set list [lsort -integer -index 0 [lsort -integer -index 0 $list]]
			set checka($chrstr) $list
		}
	}
	unset -nocomplain cachea
	if {$renamefile ne ""} {
		catch {close $f}
		set f [gzopen $renamefile]
		set header [tsv_open $f]
		set poss [list_cor $header {novelgene knowngeneid knowngene}]
		if {-1 in $poss} {error "file $renamefile should have the following fields: novelgene knowngeneid knowngene"}
		while {[gets $f line] != -1} {
			foreach {novelgene knowngeneid knowngene} [list_sub [split $line \t] $poss] break
			set cachea($novelgene) [list $knowngeneid $knowngene]
		}
	}

	# process file
	catch {sclose $f} ; catch {sclose $o}
	if {$file eq ""} {
		set f stdin
	} else {
		set f [gzopen $file]
	}
	if {$outfile eq ""} {
		set o stdout
	} else {
		set outtemp $outfile.temp[gzext $outfile]
		set o [wgzopen $outtemp]
	}
	set header [tsv_open $f]
	set poss [list_sub [tsv_basicfields $header 14 0] {0 1 2 6 12 13 7 8}]
	puts $o [join $header \t]
	set genepos [lindex $poss 4]
	set geneidpos [lindex $poss 5]

	while 1 {
		if {[gets $f line] == -1} break
		set split [split $line \t]
		if {![llength $split]} continue
		foreach {chr begin end strand gene geneid exonStarts exonEnds} [list_sub $split $poss] break
		# putsvars chr begin end strand gene geneid
		if {[regexp ^novelg_ $geneid]} {
			if {![info exists cachea($geneid)]} {
				if {![info exists checka($chr,$strand)] || $renamefile ne ""} {
					set cachea($geneid) {}
				} else {
					set list $checka($chr,$strand)
					set overlaps {}
					list_foreach {tbegin tend tgeneid tgene texonStarts texonEnds} $list {
						if {$tbegin >= $end} break
						if {$tend < $begin} continue
						# putsvars tbegin tend tgeneid tgene texonStarts texonEnds
						set obegin [max $begin $tbegin]
						set oend [min $end $tend]
						set overlap [expr {$oend - $obegin}]
						set ipct [expr {100.0*$overlap/($end-$begin)}]
						if {$ipct < 5} continue
						set opct [expr {100.0*$overlap/($tend-$tbegin)}]
						lappend overlaps [list $tgeneid $tgene $overlap $ipct $opct $obegin $oend $texonStarts $texonEnds]
					}
					# join $overlaps \n
					if {[llength $overlaps] > 1} {
						set overlaps [checkexonoverlaps $overlaps $exonStarts $exonEnds]
						if {[llength $overlaps] > 1} {
							set overlaps [lsort -index 4 -real -decreasing $overlaps]
						}
						set cachea($geneid) [lrange [lindex $overlaps 0] 0 1]
					} elseif {[llength $overlaps] == 1} {
						set cachea($geneid) [lrange [lindex $overlaps 0] 0 1]
						puts stderr "$geneid ($chr:$begin-$end) -> $cachea($geneid)"
					} else {
						set cachea($geneid) {}
					}
				}
			}
			if {[llength $cachea($geneid)]} {
				lset split $geneidpos [lindex $cachea($geneid) 0]
				lset split $genepos [lindex $cachea($geneid) 1]
				set line [join $split \t]
			}
		}
		puts $o $line
	}

	if {$o ne "stdout"} {
		close $o
		file rename -force -- $outfile.temp $outfile
	}
	if {$f ne "stdout"} {gzclose $f}

	set o [open $outfile.renames w]
	puts $o novelgene\tknowngeneid\tknowngene
	set num 0
	foreach name [bsort [array names cachea]] {
		if {$cachea($name) eq ""} continue
		incr num
		puts "$name -> $cachea($name)"
		foreach {geneid gene} $cachea($name) break
		puts $o $name\t$geneid\t$gene
	}
	puts "$num transcripts reassigned"
	close $o

}
