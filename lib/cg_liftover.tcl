proc cg_liftover {args} {
	set pos 0
	set dbdir {}
	set regionfile {}
	set correctvars 0
	set split 1
	foreach {key value} $args {
		switch -- $key {
			-dbdir {
				set dbdir $value
			}
			-regionfile - -r {
				set regionfile $value
			}
			-correctvars - -c {
				set correctvars $value
			}
			-split - -s {
				set split $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 3)} {
		errorformat liftover
		exit 1
	}
	foreach {varfile resultfile liftoverfile} $args break
	if {[file exists $resultfile]} {
		error "file $resultfile already exists, format is (now): cg liftover varfile resultfile liftoverfile"
	}

	set liftoverfile [liftoverfile $liftoverfile]
	if {[file isdir $varfile]} {
		cg_liftoversample {*}$args
	}
	set liftoverchangefile [file root $liftoverfile].refchanges.tsv
	if {[file exists $liftoverchangefile]} {
		# todo
	} else {
		# todo
	}
	set unmappedfile $resultfile.unmapped
	catch {gzclose $f} ; catch {gzclose $fl} ; catch {close $o} ; catch {close $ou}
	# open liftover file
	set fl [gzopen $liftoverfile]
	set lheader [tsv_open $fl comment]
	set lposs [list_cor $lheader {chromosome begin end strand destchromosome destbegin destend deststrand}]
	set lline [list_sub [split [gets $fl] \t] $lposs]
	set fromloc [lrange $lline 0 2]
	foreach {srcchromosome srcbegin srcend srcstrand destchromosome destbegin destend deststrand} $lline break
	if {-1 in $lposs} {exiterror "error in liftoverfile ($liftoverfile): wrong header"}
	# open varfile
	set f [gzopen $varfile]
	set header [tsv_open $f comment]
	set poss [tsv_basicfields $header 3]
	set strandpos [lsearch $header strand]
	if {$strandpos != -1} {
		lappend poss $strandpos
	}
	set strand {}
	set newheader [list_union [list_sub $header $poss] $header]
	lappend newheader beforeliftover
	# open resultfile
	set o [open $resultfile.temp w]
	puts $o "# liftover from $varfile"
	puts $o "# using $liftoverfile"
	puts $o [join $newheader \t]
	# open unmappedfile
	set ou [open $unmappedfile.temp w]
	puts $ou "# unmapped by liftover from $varfile"
	puts $ou "# using $liftoverfile"
	puts $ou [join $header \t]
	set ldone 0
	while 1 {
		if {[gets $f oline] == -1} break
		set line [split $oline \t]
		set loc [list_sub $line $poss]
		foreach {chromosome begin end strand} $loc break
		set restline [list_sub $line -exclude $poss]
		set before ${chromosome}-${begin}-${end}
		if {$strand ne ""} {append before -$strand}
		lappend restline $before
		set restline [join $restline \t]
		while 1 {
			set comp [reg_compare $fromloc $loc]
			if {$comp >= 0} {
				break
			}
			if {$ldone} break
			if {[gets $fl lline] == -1} {
				set ldone 1
				break
			}
			set lline [list_sub [split $lline \t] $lposs]
			foreach {srcchromosome srcbegin srcend srcstrand destchromosome destbegin destend deststrand} $lline break
			set fromloc [lrange $lline 0 2]
		}
		if {$comp == 0 && $begin >= $srcbegin && $end <= $srcend} {
			if {$srcstrand eq $deststrand} {
				set ustrand $strand
				set ubegin [expr {$begin + $destbegin - $srcbegin}]
				set uend [expr {$end + $destbegin - $srcbegin}]
			} else {
				if {$strand eq "+"} {set ustrand "-"} else {set ustrand "+"}
				set uend [expr {$destend - $begin + $srcbegin}]
				set ubegin [expr {$destend - $end + $srcbegin}]
			}
			if {$strandpos != -1} {
				puts $o $destchromosome\t$ubegin\t$uend\t$ustrand\t$restline
			} else {
				puts $o $destchromosome\t$ubegin\t$uend\t$restline
			}
		} else {
			puts $ou $oline
		}
	}
	catch {gzclose $f} ; catch {gzclose $fl} ; catch {close $o} ; catch {close $ou}
	#
	# sort result
	cg select -s {chromosome begin end beforeliftover} $resultfile.temp $resultfile.temp2
	file rename -force $resultfile.temp2 $resultfile
	file delete -force $resultfile.temp
	#
	# rename result, cleanup
	#
	file rename -force $unmappedfile.temp $unmappedfile
}
