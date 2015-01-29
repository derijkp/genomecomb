proc cg_liftover {args} {
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
		errorformat liftover
		exit 1
	}
	foreach {varfile liftoverfile resultfile} $args break
	if {[file exists $resultfile]} {
		error "file $resultfile already exists"
	}
	set unmappedfile $resultfile.unmapped
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set line [split [gets $f] \t]
	set poss [tsv_basicfields $header 3]
	if {![regexp ^chr [lindex $line [lindex $poss 0]]]} {set addchr 1} else {set addchr 0}
	gzclose $f
	
	# if {$poss ne "0 1 2 3 4 5"} {error "Rearranged header not supported yet, start header should be: chromosome begin end type ref alt"}
	#
	# make input file ($resultfile.temp) for liftover
	#
	set fields [list_sub $header $poss]
	set id {}
	foreach field $fields {
		lappend id \$\{$field\}
	}
	foreach {chrfield bfield efield} $fields break
	set fields {}
	if {$addchr} {
		lappend fields "chromosome=\"chr\$$chrfield\""
	} else {
		lappend fields "chromosome=\"\$$chrfield\""
	}
	lappend fields "begin=if(\$$bfield == \$$efield,\$$bfield - 1,\$$bfield)"
	lappend fields "end=if(\$$bfield == \$$efield,\$$efield + 1,\$$efield)"
	lappend fields id=\"[join $id -]\"
	cg select -f $fields -sh $resultfile.temph $varfile $resultfile.temp
	#
	# do liftover -> $resultfile.temp2
	#
	# set dir [file dir [exec which liftOver]]
	if {[catch {exec liftOver -bedPlus=3 -tab $resultfile.temp $liftoverfile $resultfile.temp2 $unmappedfile.temp} errmsg]} {
		puts "$errmsg"
	}
	#
	# add original data back to liftovered file -> $resultfile.temp3
	# we cannot just let liftover do it, as it only takes max 12 columns with it
	#
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set fl [open $resultfile.temp2]
	#set temp [tsv_open $fl]
	set o [open $resultfile.temp3 w]
	set header [list_union [list_sub $header $poss] $header]
	lappend header beforeliftover
	puts $o "# liftover from $varfile"
	puts $o "# using $liftoverfile"
	puts $o [join $header \t]
	while {![eof $fl]} {
		set lline [split [gets $fl] \t]
		if {![llength $lline]} continue
		foreach {chr begin end lname} $lline break
		foreach {lchr lbegin lend} [split $lname -] break
		while 1 {
			set line [split [gets $f] \t]
			set loc [list_sub $line $poss]
			set name [join $loc -]
			if {$name eq $lname} break
			if {[eof $f]} {error "$lname not found"}
		}
		if {$lbegin == $lend} {
			incr begin ; incr end -1
		}
		puts $o $chr\t$begin\t$end\t[join [list_sub $line -exclude $poss] \t]\t$lname
	}
	close $o
	close $fl
	gzclose $f
	#
	# sort result -> $resultfile.temp4
	#
	cg select -s - $resultfile.temp3 $resultfile.temp4
	#
	# rename result, cleanup
	#
	file rename -force $resultfile.temp4 $resultfile
	file delete $resultfile.temph $resultfile.temp $resultfile.temp2 $resultfile.temp3
	#
	# fo unmapped: add original data back
	#
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set fl [open $unmappedfile.temp]
	#set temp [tsv_open $fl]
	set o [open $unmappedfile.temp3 w]
	set header [list_union [list_sub $header $poss] $header]
	puts $o "# unmapped by liftover from $varfile"
	puts $o "# using $liftoverfile"
	puts $o [join $header \t]
	while {![eof $fl]} {
		set lline [gets $fl]
		if {[string index $lline 0] eq "\#"} continue
		set lline [split $lline \t]
		if {![llength $lline]} continue
		foreach {chr begin end lname} $lline break
		foreach {lchr lbegin lend} [split $lname -] break
		while 1 {
			set line [split [gets $f] \t]
			set loc [list_sub $line $poss]
			set name [join $loc -]
			if {$name eq $lname} break
			if {[eof $f]} {error "$lname not found"}
		}
		if {$lbegin == $lend} {
			incr begin ; incr end -1
		}
		puts $o $chr\t$begin\t$end\t[join [list_sub $line -exclude $poss] \t]
	}
	close $o
	close $fl
	gzclose $f
	#
	# sort result -> $unmappedfile.temp4
	#
	cg select -s - $unmappedfile.temp3 $unmappedfile.temp4
	#
	# rename result, cleanup
	#
	file rename -force $unmappedfile.temp4 $unmappedfile
	file delete $unmappedfile.temp2 $unmappedfile.temp3
}
