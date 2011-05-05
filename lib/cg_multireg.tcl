#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

if 0 {
	package require Tclx
	signal -restart error SIGINT
	lappend auto_path /home/peter/dev/completegenomics/lib
	package require Extral
	cd /complgen/multireg
	set compar_file multireg.tsv
	set file ../GS102/sreg-GS102.tsv
	set file ../GS103/sreg-GS103.tsv

	~/dev/completegenomics/lib/cg_multireg.tcl multireg.tsv ../GS102/sreg-GS102.tsv
	~/dev/completegenomics/lib/cg_multireg.tcl testmultireg.tsv sreg-test.tsv
}


proc multireg_next1 {f1 poss1} {
	set cur1 {}
	while {![eof $f1]} {
		set cur1 [split [gets $f1] \t]
		if {![llength $cur1]} continue
		set part1 [lrange $cur1 3 end]
		set comp1 [list_sub $cur1 $poss1]
		foreach {chr1 start1 end1} $comp1 break
		if {[isint $start1]} break
	}
	if {![llength $cur1] || ![isint $start1]} {return {{} 1000 1000 -1 -1}}
	set nchr1 [chr2num $chr1]
	list $part1 $chr1 $nchr1 $start1 $end1
}

proc multireg_next2 {f2 poss2} {
	set start2 {}
	while {![eof $f2]} {
		set comp2 [get_region $f2 $poss2]
		foreach {chr2 start2 end2} $comp2 break
		if {[isint $start2]} break
	}
	if {![isint $start2]} {return {1000 1000 -1 -1}}
	set nchr2 [chr2num $chr2]
	list $chr2 $nchr2 $start2 $end2
}

proc multireg {compar_file file} {
	global cache comparposs1 mergeposs1 comparposs2 mergeposs2 dummy1 dummy2 restposs1 restposs2

	set name [file root [file tail [gzfile $file]]]
	catch {close $f1}; catch {close $f2}; catch {close $o}
	set f2 [gzopen $file]
	set poss2 [open_region $f2 h2]
	set num 0
	if {![file exists $compar_file]} {
		close $f2
		set h2base [list_sub $h2 $poss2]
		cg select \
			-f "chromosome=\$[lindex $h2base 0] begin=\$[lindex $h2base 1] end=\$[lindex $h2base 2] $name=1" \
			$file $compar_file
		return
	}
	set f1 [open $compar_file]
	set poss1 [open_region $f1 h1]
	if {[inlist $h1 $name]} {
		error "$name already present in $compar_file"
	}
	set o [open $compar_file.temp w]
	puts $o [join $h1 \t]\t$name
	set dummy1 [list_fill [expr {[llength $h1]-3}] 0]
	close $o
	# putslog "multireg $compar_file $poss1 $dummy1 $file $poss2 >> $compar_file.temp"
	exec multireg $compar_file {*}$poss1 [join $dummy1 \t] $file {*}$poss2 >> $compar_file.temp 2>@stderr
	catch {file rename -force $compar_file $compar_file.old}
	file rename $compar_file.temp $compar_file	
}

proc cg_multireg {args} {
	if {([llength $args] < 1)} {
		errorformat multireg
		exit 1
	}
	foreach {compar_file} $args break
	if {[file exists $compar_file]} {
		set f [gzopen $compar_file]
		set header [tsv_open $f]
		close $f
	} else {
		set header {chromosome begin end}
	}
	set files [lrange $args 1 end]
	foreach file $files {
		set name [file root [file tail [gzfile $file]]]
		if {[inlist $header $name]} {
			putslog "*** Skipping $file: $name already in $compar_file ***"
			continue
		}
		putslog "Adding $file to $compar_file"
		multireg $compar_file $file
	}
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base [file tail [info script]]
	cg_multireg {*}$argv
}
