proc cg_exportplink {args} {
	set query {}
	set pos 0
	set codegeno 0
	foreach {key value} $args {
		switch -- $key {
			-q {
				set query $value
			}
			-c - --codegeno {
				set codegeno $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] != 2} {
		errorformat cg_exportplink
		exit 1
	}
	foreach {varfile resultfile} $args break
	catch {close $f} ; catch {close $o}
	if {$query ne ""} {
		set f [open "|[list cg select -q $query $varfile]"]
	} else {
		set f [gzopen $varfile]
	}
	set header [tsv_open $f]
	set alleleposs [list_find -glob $header alleleSeq*]
	set temp {}
	foreach field [list_sub $header $alleleposs] p $alleleposs {
		regexp {^([^-]+)-(.*)$} $field unused an name
		lappend temp [list $name $an $p]
	}
	set temp [ssort -natural $temp]
	set o [open $resultfile.tfam.pre w]
	foreach {name ignore} [list_subindex $temp 0] {
		puts $o [join [list fam $name 0 0 0 -9] \t]
	}
	close $o
	set aposs [list_subindex $temp 2]
	set poss [tsv_basicfields $header 6]
	set o [open $resultfile.tped w]
#	array set trans {? 0 - 0 N 0 {} -}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line] < 4} continue
		foreach {chr b e type ref alts} [list_sub $line $poss] break
		if {![isint $b]} {
			puts stderr "skipping var: error in line: $line"
			continue
		}
		set chr [chr_clip $chr]
		if {$chr eq "M"} {set chr MT}
		set alts [split $alts ,]
		if {$codegeno} {
			set refcode 1
		} else {
			set refcode $ref
		}
		if {[llength $alts] > 1} {
			puts stderr "Warning: more than 2 alleles for $chr-$b-$e-$type: splitting"
		}
		foreach alt $alts {
			if {$codegeno} {
				set altcode 2
			} else {
				set altcode $alt
			}
			set result [list $chr $chr-$b-$e-$type-$alt [format %.4f [expr {$b/1000000.0}]] $b]
			foreach {gt1 gt2} [list_sub $line $aposs] {
				if {$gt1 eq $ref} {
					set gt1 $refcode
				} elseif {$gt1 eq $alt} {
					set gt1 $altcode
				} else {
					set gt1 0
				}
				if {$gt2 eq $ref} {
					set gt2 $refcode
				} elseif {$gt2 eq $alt} {
					set gt2 $altcode
				} else {
					set gt2 0
				}
				if {$gt1 eq "0"} {set gt2 0}
				if {$gt2 eq "0"} {set gt1 0}
				lappend result $gt1 $gt2
			}
			puts $o [join $result \t]
		}
	}
	close $o
	catch {close $f}
}

