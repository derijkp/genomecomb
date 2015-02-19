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
	if {[file exists $resultfile]} {
		error "file $resultfile already exists"
	}
	set unmappedfile $resultfile.unmapped
	set tempfile [scratchfile]
	set unmappedtempfile [scratchfile]
	set tempfile2 [scratchfile]
	set tempfile3 [scratchfile]
	set tempfile4 [scratchfile]

	catch {close $f} ; catch {close $fl} ; catch {close $o}
	set f [gzopen $file]
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 3]
	set line [split [gets $f] \t]
	if {![regexp ^chr [lindex $line [lindex $poss 0]]]} {set addchr 1} else {set addchr 0}
	#
	# make input file ($resultfile.temp) for liftover
	# split region into snps, that will be lift over separately
	#
	set chrpos [lindex $poss 0]
	set o [open $tempfile w]
	while {![eof $f]} {
		set part [list_sub $line $poss]
		foreach {chr begin end} $part break
		if {$addchr} {set chr chr$chr}
		set ori ${chr}-$begin-$end
		while {$begin < $end} {
			puts $o $chr\t$begin\t[incr begin]\t$ori
		}
		set line [split [gets $f] \t]
		if {![llength $line] && [eof $f]} break
	}
	close $o
	close $f
	#
	# do liftover -> $resultfile.temp2
	#
	# set dir [file dir [exec which liftOver]]
	if {[catch {exec liftOver -bedPlus=3 -tab $tempfile $liftoverfile $tempfile2 $unmappedtempfile} errmsg]} {
		puts "$errmsg"
	}
	#
	# rejoin regions that were separated as far as possible
	# add original data back to liftovered file -> $tempfile3
	# we cannot just let lifover do it, as it only takes max 12 columns with it
	#
	set f [gzopen $file]
	set header [tsv_open $f comment]
	set fl [open $tempfile2]
	set o [open $tempfile3 w]
	set header [list_union [list_sub $header $poss] $header]
	lappend header beforeliftover
	puts $o [join $header \t]
	set prevname {} ; set prevchr {} ; set prevend {} ; set  prevbegin {}
	set rising {}
	while {![eof $fl]} {
		set lline [split [gets $fl] \t]
		if {![llength $lline]} continue
		foreach {chr begin end lname} $lline break
		set break 0
		if {$begin eq $prevend} {
			if {$rising eq ""} {
				set rising 1
			} elseif {$rising eq "0"} {
				set rising {}
				set break 1
			}
		} elseif {$end eq $prevbegin} {
			if {$rising eq ""} {
				set rising 0
			} elseif {$rising eq "1"} {
				set rising {}
				set break 1
			}
		} else {
			set break 1
		}
		
		if {$chr ne $prevchr || $break || $lname ne $prevname} {
			if {$prevname ne ""} {
				puts $o $prevchr\t$prevbegin\t$prevend\t[join [list_sub $line -exclude $poss] \t]\t$prevname
			}
			if {$lname ne $prevname} {
				while 1 {
					set line [split [gets $f] \t]
					if {$addchr} {
						set name chr[join [list_sub $line $poss] -]
					} else {
						set name [join [list_sub $line $poss] -]
					}
					if {$name eq $lname} break
					if {[eof $f]} {error "$lname not found"}
				}
			}
			set prevchr $chr
			set prevbegin $begin
			set prevend $end
			set prevname $lname
			set rising {}
		} else {
			if {$rising eq "0"} {
				set prevbegin $begin
			} else {
				set prevend $end
			}
		}
	}
	close $o
	close $fl
	gzclose $f
	#
	# sort result -> $tempfile4
	#
	cg select -s - $tempfile3 $tempfile4
	# collapse overlapping regions
	set o [open $resultfile.temp w]
	puts $o "# liftregion from $file"
	puts $o "# using $liftoverfile"
	close $o
	cg regcollapse $tempfile4 >> $resultfile.temp
	#
	# rename result, cleanup
	#
	file rename -force $resultfile.temp $resultfile
	file delete $tempfile $tempfile2 $tempfile3

	#
	# rejoin regions of unmapped file that were separated as far as possible
	#
	set f [gzopen $file]
	set header [tsv_open $f comment]
	set fl [open $unmappedtempfile]
	set o [open $tempfile3 w]
	set header [list_union [list_sub $header $poss] oriregion $header]
	puts $o [join $header \t]
	set prevname {} ; set prevchr {} ; set prevend {} ; set  prevbegin {}
	set rising {}
	while {![eof $fl]} {
		set lline [gets $fl]
		if {[string index $lline 0] eq "\#"} continue
		if {$lline eq ""} continue
		set lline [split $lline \t]
		foreach {chr begin end lname} $lline break
		if {$chr ne $prevchr || $lname ne $prevname || $begin ne $prevend} {
			if {$prevname ne ""} {
				puts $o $prevchr\t$prevbegin\t$prevend\t$prevname\t[join [list_sub $line -exclude $poss] \t]\t
			}
			if {$lname ne $prevname} {
				while 1 {
					set line [split [gets $f] \t]
					if {$addchr} {
						set name chr[join [list_sub $line $poss] -]
					} else {
						set name [join [list_sub $line $poss] -]
					}
					if {$name eq $lname} break
					if {[eof $f]} {error "$lname not found"}
				}
			}
			set prevchr $chr
			set prevbegin $begin
			set prevend $end
			set prevname $lname
		} else {
			set prevend $end
		}
	}
	close $o
	close $fl
	gzclose $f
	#
	# sort result -> $tempfile4
	#
	cg select -s - $tempfile3 $tempfile4
	# collapse overlapping regions
	set o [open $unmappedfile.temp w]
	puts $o "# unmapped liftregion from $file"
	puts $o "# using $liftoverfile"
	close $o
	cg regcollapse $tempfile4 >> $unmappedfile.temp
	#
	# rename result, cleanup
	#
	file rename -force $unmappedfile.temp $unmappedfile
	file delete $tempfile $unmappedtempfile $tempfile2 $tempfile3
}
