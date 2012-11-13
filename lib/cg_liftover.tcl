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
	if {![regexp ^chr [lindex $line 0]]} {set addchr 1} else {set addchr 0}
	gzclose $f
	set poss [tsv_basicfields $header 3]
	# if {$poss ne "0 1 2 3 4 5"} {error "Rearranged header not supported yet, start header should be: chromosome begin end type ref alt"}
	#
	# make input file ($resultfile.temp) for liftover
	#
	set fields [list_sub $header $poss]
	set id {}
	foreach field $fields {
		lappend id \$\{$field\}
	}
	if {$addchr} {
		set fields [lreplace $fields 0 0 "chromosome=\"chr\$chromosome\""]
	}
	lappend fields id=\"[join $id -]\"
	cg select -f $fields -sh $resultfile.temph $varfile $resultfile.temp
	#
	# do liftover -> $resultfile.temp2
	#
	# set dir [file dir [exec which liftOver]]
	if {[catch {exec liftOver -bedPlus=3 -tab $resultfile.temp $liftoverfile $resultfile.temp2 $unmappedfile} errmsg]} {
		puts "$errmsg"
	}
	#
	# add original data back to liftovered file -> $resultfile.temp3
	# we cannot just let lifover do it, as it only takes max 12 columns with it
	#
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set fl [open $resultfile.temp2]
	set temp [tsv_open $fl]
	set o [open $resultfile.temp3 w]
	lappend header beforeliftover
	puts $o "# liftover from $varfile"
	puts $o "# using $liftoverfile"
	puts $o [join $header \t]
	while {![eof $fl]} {
		set lline [split [gets $fl] \t]
		if {![llength $lline]} continue
		set lname [lindex $lline 3]
		while 1 {
			set line [split [gets $f] \t]
			set name [join [list_sub $line $poss] -]
			if {$name eq $lname} break
			if {[eof $f]} {error "$lname not found"}
		}
		set rline [lrange $lline 0 2]
		lappend rline {*}[list_sub $line -exclude $poss] [join [lrange [split $lname -] 0 2] -]
		puts $o [join $rline \t]
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
	file rename $resultfile.temp4 $resultfile
	file delete $resultfile.temph $resultfile.temp $resultfile.temp2 $resultfile.temp3
}
