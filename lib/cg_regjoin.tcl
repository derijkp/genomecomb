proc regjoin_fields {regfile1 fields} {
	catch {close $f} ; catch {close $o}
	if {$regfile1 ne "stdin"} {
		if {[catch {gzopen $regfile1} f]} {
			error "Could not open $regfile1"
		}
	} else {
		set f stdin
	}
	set header [tsv_open $f]
	set cor [tsv_basicfields $header 3]
	foreach {chrpos startpos endpos} $cor break
	if {$fields eq "*"} {
		set fields [list_sub $header -exclude $cor]
	} else {
		set fields [list_common $fields $header]
	}
	lappend cor {*}[list_cor $header $fields]
	set nh [list_sub $header $cor]
	puts [join $nh \t]
	set pline [split [gets $f] "\t"]
	set pline [list_sub $pline $cor]
	foreach {pchr pstart pend} $pline break
	while {![eof $f]} {
		if {[gets $f line] == -1} break
		set line [split $line "\t"]
		set line [list_sub $line $cor]
		foreach {chr start end} $line break
		set join 1
		if {$chr ne $pchr} {
			set join 0
		} elseif {$start > $pend} {
			set join 0
		} elseif {[lrange $line 3 end] ne [lrange $pline 3 end]} {
			set join 0
		}
		if {!$join} {
			puts [join $pline \t]
			set pline $line
			foreach {pchr pstart pend} $pline break
		} else {
			lset pline 2 $end
			set pend $end
		}
	}
	puts [join $pline \t]
	close $f
}

proc regjoin {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	if {$regfile2 ne ""} {
		set f2 [gzopen $regfile2]
		set poss2 [open_region $f2]
		gzclose $f2
	} else {
		set poss2 {0 0 0}
	}
	if {$regfile1 ne ""} {
		set f1 [gzopen $regfile1]
		set poss1 [open_region $f1]
		gzclose $f1
		set temp1 [gztemp $regfile1]
		set temp2 [gztemp $regfile2]
		# puts [list reg_join $regfile1 {*}$poss1 $regfile2 {*}$poss2]
		exec reg_join $temp1 {*}$poss1 $temp2 {*}$poss2 >@ stdout 2>@ stderr
#		gzrmtemp $temp1 $temp2
	} else {
		set poss1 [open_region stdin]
		set temp2 [gztemp $regfile2]
		set o [open "| [list reg_join - {*}$poss1 $temp2 {*}$poss2] >@ stdout 2>@ stderr" w]
		fcopy stdin $o
		close $o
		gzrmtemp $temp2
	}
}

proc cg_regjoin {args} {
	foreach {region_file1 region_file2} {{} {}} break
	set fields {}
	cg_options regjoin args {
		-fields {
			set fields $value
		}
	} {region_file1 region_file2} 0 2
	if {$fields ne ""} {
		if {$region_file2 ne ""} {error "regjoin with -fields only supported for one file"}
		regjoin_fields $region_file1 $fields
	} else {
		regjoin $region_file1 $region_file2
	}
}
