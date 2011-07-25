proc cg_exportplink {args} {
	if {[llength $args] == 0} {
		errorformat cg_export_plink
		exit 1
	}

	set query {}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-q {
				set query $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	foreach {varfile resultfile} $args break
	if {$query ne ""} {
		set f [open "|[list cg select -q $query $varfile]"]
	} else {
		set f [gzopen $varfile]
	}
	set header [tsv_open $f]
	set alleleposs [list_find -glob $header alleleSeq*]
	set temp {}
	foreach field [list_sub $header $alleleposs] p $alleleposs {
		foreach {an name} [split $field -] break
		lappend temp [list $name $an $p]
	}
	set temp [lsort -dict $temp]
	set o [open $resultfile.tfam.pre w]
	foreach {name ignore} [list_subindex $temp 0] {
		puts $o [join [list fam $name 0 0 0 -9] \t]
	}
	close $o
	set aposs [list_subindex $temp 2]
	set poss [tsv_basicfields $header 4]
	set o [open $resultfile.tped w]
	array set trans {? 0 - 0 N 0 {} -}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line] < 4} continue
		foreach {chr b e type} [list_sub $line $poss] break
		if {![isint $b]} {
			puts stderr "skipping var: error in line: $line"
			continue
		}
		regsub ^chr $chr {} chr
		if {$chr eq "M"} {set chr MT}
		set result [list $chr $chr-$b-$e-$type [format %.4f [expr {$b/1000000.0}]] $b]
		set alleles {}
		unset -nocomplain a
		foreach {gt} [list_sub $line $aposs] {
			set gt [get trans($gt) $gt]
			lappend alleles $gt
			incr a($gt)
		}
		unset -nocomplain a(0)
		set names [array names a]
		if {[llength $names] > 2} {
			set temp {}
			foreach allele $names {
				if {$type eq "del" && ($allele eq "-")} continue
				lappend temp [list $a($allele) $allele]
			}
			set temp [lsort -index 0 -integer -decreasing $temp]
			if {$type eq "del"} {set s 1} else {set s 2}
			set remove {}
			foreach e [list_subindex [lrange $temp $s end] 1] {
				lappend remove $e 0
			}
			set alleles [list_change $alleles $remove]
			puts stderr "Warning: more than 2 alleles for $result: converted $remove"
			continue
		}
		foreach {gt1 gt2} $alleles {
			if {$gt1 eq "0"} {set gt2 0}
			if {$gt2 eq "0"} {set gt1 0}
			lappend result $gt1 $gt2
		}
		puts $o [join $result \t]
	}

	close $o
	catch {close $f}
}

