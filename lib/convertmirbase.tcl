proc convertmirbase_findloop {resultline mirnas seq structs} {
	foreach {chr begin end strand name} $resultline break
	if {[dict exists $structs $name]} {
		foreach {loopbegin looplen} [dict get $structs $name] break
	} else {
		set rnafold [split [exec [rnafold] --noPS --MEA -p -d2 << $seq] \n]
		set struct [lindex [lindex $rnafold 1] 0]
		foreach {loopbegin loopend} [findloop $struct] break
	}
	if {$strand eq "+"} {
		set loopbegin [expr {$begin + $loopbegin}]
		set loopend [expr {$loopbegin + $looplen}]
	} else {
		set loopend [expr {$end - $loopbegin}]
		set loopbegin [expr {$loopend - $looplen}]
	}
	set mirnas [lsort -index 0 -integer $mirnas]
	if {[llength $mirnas] > 1} {
		foreach {b1 e1 b2 e2} [list_concat $mirnas] break
		if {$e1 > $b2} {
			error "overlapping mirnas [join $mirnas ,] for $resultline"
		}
		if {$loopbegin >= $b2 || $loopend <= $e1} {
			log "[lindex $resultline 4]: loop outside of mirnas, shifted"
			set loopbegin $b2
			set loopend $e2
		}
#		if {$loopbegin < $e1} {set loopbegin $e1}
#		if {$loopend > $b2} {set loopend $b2}
		lset resultline 5 $b1
		lset resultline 6 $e1
		lset resultline 9 $b2
		lset resultline 10 $e2
	} elseif {[llength $mirnas] == 1} {
		list_foreach {b e} $mirnas break
		if {$e > $loopend && $b < $loopbegin} {
			putslog "[lindex $resultline 4]: loop overlaps mirna, shifted"
			if {[expr {$b-$loopbegin}] < [expr {[lindex $resultline 2]-$e}]} {
				set loopbegin $e ; set loopend [expr {$e+4}]
			} else {
				set loopbegin [expr {$b-4}] ; set loopend $b
			}
		}
		if {$e <= $loopend} {
			lset resultline 5 $b
			lset resultline 6 $e
#			if {$loopbegin < $e} {set loopbegin $e}
		} else {
			lset resultline 9 $b
			lset resultline 10 $e
#			if {$loopend > $b} {set loopend $b}
		}
	}
	lset resultline 7 $loopbegin
	lset resultline 8 $loopend
	return $resultline
}

proc convertmirbase {gff3file resultfile genomefile {structfile {}}} {
# putsvars gff3file resultfile genomefile structfile
	if {[file exists $resultfile]} {error "File $resultfile exists"}
	set gf [genome_open $genomefile]
	set structs [dict create]
	if {$structfile ne ""} {
		set f [open $structfile]
		while {![eof $f]} {
			set line [gets $f]
			if {[regexp {^>([^ ]+)} $line temp name]} {
				set line [gets $f]
				set bulges1 [gets $f]
				set arm1 [gets $f]
				set mid [gets $f]
				set arm2 [gets $f]
				set bulges2 [gets $f]
				set line [gets $f]
				set helixlast [string last | $mid]
				set arm1bases [string range $bulges1 0 $helixlast][string range $arm1 0 $helixlast]
				set loopbegin [regsub -all {[^ |-]} $arm1bases {} temp]
				incr helixlast
				set looplen 0
				foreach string [list $bulges1 $arm1 $mid $arm2 $bulges2] {
					set temp [string range $string $helixlast end]
					incr looplen [regsub -all {[^ |-]} $temp {} temp]
				}
				dict set structs $name [list $loopbegin $looplen]
			}
		}
		close $f
	}
	catch {close $f}
	set f [open $gff3file]
	set id {}
	set result {}
	set resultline {}
	set mirnas {}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] eq "#"} continue
		set line [split $line \t]
		foreach {chr source type begin end score strand phase attr} $line break
		incr begin -1
		if {![llength $line] || $type eq "miRNA_primary_transcript"} {
			if {[llength $resultline]} {
				# add previuous to results
				if {$seq ne ""} {
					set resultline [convertmirbase_findloop $resultline $mirnas $seq $structs]
					lappend result $resultline
				}
				set mirnas {}
			}
			if {![llength $line]} continue
			if {![regexp {Name=([^;\n]+)} $attr temp name]} {
				error "Could not get name from $line"
			}
			if {![regexp {ID=([^;\n]+)} $attr temp id]} {
				error "Could not get id from $line"
			}
			if {[catch {
				set seq [string toupper [genome_get $gf $chr $begin $end]]
			} msg]} {
				puts stdout "skipping $name: Could not find $chr:${begin}-$end in genome"
				set seq {}
			}
			set resultline [list $chr $begin $end $strand $name {} {} {} {} {} {}]
		} elseif {$type eq "miRNA"} {
			if {![regexp {Derives_from=([^;\n]+)} $attr temp fromid]} {
				error "Could not get Derives_from from $line"
			}
			if {$fromid ne $id} {
				if {![string match ${fromid}_* $id]} {
					error "miRNA ($line) does not derive from $id: $resultline"
				} else {
					puts "warning: approx match for miRNA ($line) derive from $id: $resultline"
				}
			}
			lappend mirnas [list $begin $end]
		}
	}
	close $f
	set result [lsort -integer -index 1 $result]
	set result [lsort -dict -index 0 $result]
	set o [open $resultfile.temp w]
	puts $o [join {chromosome begin end strand name mature1start mature1end loopstart loopend mature2start mature2end} \t]
	foreach resultline $result {
		puts $o [join $resultline \t]
	}
	close $o
	if {[file extension $resultfile] eq ".lz4"} {
		cg lz4 -i 1 $resultfile.temp
		file rename -force $resultfile.temp.lz4 $resultfile
		file rename -force $resultfile.temp.lz4.lz4i $resultfile.lz4i
	} else {
		file rename -force $resultfile.temp $resultfile
	}
}
