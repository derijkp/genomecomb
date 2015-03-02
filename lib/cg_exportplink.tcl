proc cg_exportplink {args} {
	set query {}
	set pos 0
	set nulllines 0
	set codegeno 0
	foreach {key value} $args {
		switch -- $key {
			-q {
				set query $value
			}
			-c - --codegeno {
				set codegeno $value
			}
			-s - --samples {
				set samples $value
			}
			-n - --nulllines {
				set nulllines $value
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
	if {![llength $header]} {
		catch {close $f} msg
		error "error querying file: $msg"
	}
	if {![info exists samples]} {
		set samples [samples $header]
	}
	set aposs {}
	set samplepresents {}
	foreach sample $samples {
		if {$sample eq ""} {
			lappend samplepresents 0
			lappend aposs -1 -1
			continue
		} else {
			lappend samplepresents 1
		}
		set pos [lsearch $header alleleSeq1-$sample]
		if {$pos == -1} {error "no genotype (alleleSeq1-$sample) found for $sample"}
		lappend aposs $pos
		set pos [lsearch $header alleleSeq2-$sample]
		if {$pos == -1} {error "no genotype (alleleSeq2-$sample) found for $sample"}
		lappend aposs $pos
	}
	set o [open $resultfile.tfam.pre w]
	foreach name $samples {
		puts $o [join [list fam $name 0 0 0 -9] \t]
	}
	close $o
	set poss [tsv_basicfields $header 6]
	set o [open $resultfile.tped w]
#	array set trans {? 0 - 0 N 0 {} -}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line] < 4} continue
		foreach {chr b e type ref alts} [list_sub $line $poss] break
		if {![isint $b]} {
			putslog "skipping var: error in line: $line"
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
		if {[llength $alts] == 0} {
			set alts {{}}
		} elseif {[llength $alts] > 1} {
			putslog "Warning: more than 2 alleles for $chr-$b-$e-$type: splitting"
		}
		foreach alt $alts {
			if {$codegeno} {
				set altcode 2
			} else {
				set altcode $alt
			}
			set result {}
			foreach {gt1 gt2} [list_sub $line $aposs] samplepresent $samplepresents {
				if {!$samplepresent} {
					lappend result 0 0
					continue
				}
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
			if {$nulllines == 0} {
				if {![llength [list_remove $result 0]]} continue
			}
			puts $o [join [list $chr $chr-$b-$e-$type-$ref-$alt [format %.6f [expr {$b/1000000.0}]] $b {*}$result] \t]
		}
	}
	close $o
	catch {close $f}
}

