proc cg_wide {args} {
	set samplefields {}
	set commonfields {}
	set file {}
	set outfile {}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-s - --samplefields {
				set samplefields $value
			}
			-f - --fields {
				set commonfields $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {[llength $args] > 2} {
		errorformat wide
	}
	foreach {file outfile} $args break
	if {$file eq ""} {
		set f stdin
		set tempfile [tempfile]
	} else {
		set f [gzopen $file]
		set tempfile {}
	}
	if {$outfile eq ""} {
		set o stdout
	} else {
		set o [open $outfile.temp w]
	}
	set header [tsv_open $f comment]
	if {$commonfields eq ""} {
		set poss [tsv_basicfields $header 6 0]
	} else {
		set commonfields [expandfields $header $commonfields]
		set poss [list_cor $header $commonfields]
	}
	if {$samplefields eq ""} {
		set samplefields [list_common $header {sample sample1 sample2 sample3 mapping varcall}]
	}
	set poss [list_remove $poss -1]
	set commonfields [list_sub $header $poss]
	set samplefields [expandfields $header $samplefields]
	set sposs [list_remove [list_cor $header $samplefields] -1]
	set sfields [list_sub $header $sposs]
	set rfields [list_sub $header -exclude [list_union $poss $sposs]]
	set rposs [list_cor $header $rfields]
	# find out samples in file (parse entirely)
	putslog "Getting samples"
	unset -nocomplain a
	set prevcur {}
	set sort 0
	if {$tempfile ne ""} {
		set tf [open $tempfile w]
		puts $tf $comment
		puts $tf [join $header \t]
		while {![eof $f]} {
			set oline [gets $f]
			set line [split $oline \t]
			if {![llength $line]} continue
			set a([list_sub $line $sposs]) 1
			set cur [list_sub $line $poss]
			if {[loc_compare $prevcur $cur] > 0} {set sort 1}
			set prevcur $cur
			puts $tf $oline
		}
		close $tf
	} else {
		while {![eof $f]} {
			set oline [gets $f]
			set line [split $oline \t]
			if {![llength $line]} continue
			set a([list_sub $line $sposs]) 1
			set cur [list_sub $line $poss]
			if {[loc_compare $prevcur $cur] > 0} {set sort 1}
			set prevcur $cur
		}
		close $f
		set tempfile $file
	}
	set samples [lsort -dict [array names a]]
	# sort ?
	if {$sort} {
		putslog "Sorting"
		set tempfile2 [tempfile]
		cg select -s $commonfields $tempfile $tempfile2
		set tempfile $tempfile2
	}
	#
	putslog "Creating wide file"
	set f [gzopen $tempfile]
	tsv_open $f
	set newheader $commonfields
	foreach sample $samples {
		foreach field $rfields {
			lappend newheader $field-[join $sample -]
		}
	}
	puts $o "#type\twide tsv"
	puts $o "#samplefields\t$sfields"
	puts -nonewline $o $comment
	puts $o [join $newheader \t]
	unset -nocomplain a
	set def \t[join [list_fill [llength $rfields] ?] \t]
	set prev {}
	set end 0
	while {1} {
		set line [split [gets $f] \t]
		if {![llength $line]} {
			if {[eof $f]} {set end 1} else continue
		}
		set cur [list_sub $line $poss]
		if {$prev eq ""} {set prev $cur}
		if {$cur ne $prev || $end} {
			puts -nonewline $o [join $prev \t]
			foreach sample $samples {
				if {[info exists a($sample)]} {
					puts -nonewline $o \t[join $a($sample) \t]
				} else {
					puts -nonewline $o $def
				}
			}
			puts $o ""
			if {$end} break
			unset -nocomplain a
			set prev $cur
		}
		set sample [list_sub $line $sposs]
		set a($sample) [list_sub $line $rposs]
	}
	if {$o ne "stdout"} {
		close $o
		file rename -force $outfile.temp $outfile
	}
	if {$f ne "stdin"} {gzclose $f}
}

if 0 {
	set tsvfile tests/data/reg4.tsv
	set tsvfile tests/data/testvars.tsv
	set o stdout
}
