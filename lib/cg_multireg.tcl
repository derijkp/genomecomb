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

	set name [file root [file tail $file]]
	regexp {[^-]+$} $name name
	catch {close $f1}; catch {close $f2}; catch {close $o}
	set f2 [open $file]
	set poss2 [open_region $f2 h2]
	set num 0
	if {![file exists $compar_file]} {
		incr num
		if {![expr {$num % 100000}]} {putslog $num}
		set o [open $compar_file w]
		puts $o [join {chromosome begin end} \t]\tsequenced-$name
		while {![eof $f2]} {
			set line [get_region $f2 $poss2]
			if {![llength $line]} continue
			puts $o [join $line \t]\t1
		}
		close $o
		return
	}
	set f1 [open $compar_file]
	set poss1 [open_region $f1 h1]
	if {[inlist $h1 sequenced-$name]} {
		error "$name already present in $compar_file"
	}
	set o [open $compar_file.temp w]
	puts $o [join $h1 \t]\tsequenced-$name
	set dummy1 [list_fill [expr {[llength $h1]-3}] 0]
	foreach {part1 chr1 nchr1 start1 end1} [multireg_next1 $f1 $poss1] break
	foreach {chr2 nchr2 start2 end2} [multireg_next2 $f2 $poss2] break
	while {![eof $f1] || ![eof $f2]} {
		incr num
		if {![expr {$num % 100000}]} {putslog $num}
		if {($nchr2 < $nchr1) || ($nchr1 == $nchr2 && $end2 <= $start1)} {
			puts $o $chr2\t$start2\t$end2\t[join $dummy1 \t]\t1
			foreach {chr2 nchr2 start2 end2} [multireg_next2 $f2 $poss2] break
		} elseif {($nchr1 < $nchr2) || ($nchr1 == $nchr2 && $end1 <= $start2)} {
			puts $o $chr1\t$start1\t$end1\t[join $part1 \t]\t0
			foreach {part1 chr1 nchr1 start1 end1} [multireg_next1 $f1 $poss1] break
		} else {
			if {$start1 < $start2} {
				puts $o $chr1\t$start1\t$start2\t[join $part1 \t]\t0
				set start1 $start2
			}
			if {$start2 < $start1} {
				puts $o $chr1\t$start2\t$start1\t[join $dummy1 \t]\t1
				set start2 $start1
			}
			if {$end2 < $end1} {
				puts $o $chr1\t$start1\t$end2\t[join $part1 \t]\t1
				set start1 $end2
				foreach {chr2 nchr2 start2 end2} [multireg_next2 $f2 $poss2] break
			} elseif {$end1 < $end2} {
				puts $o $chr1\t$start1\t$end1\t[join $part1 \t]\t1
				set start2 $end1
				foreach {part1 chr1 nchr1 start1 end1} [multireg_next1 $f1 $poss1] break
			} else {
				puts $o $chr1\t$start1\t$end1\t[join $part1 \t]\t1
				foreach {chr2 nchr2 start2 end2} [multireg_next2 $f2 $poss2] break
				foreach {part1 chr1 nchr1 start1 end1} [multireg_next1 $f1 $poss1] break
			}
		}
	}

	close $f1; close $f2; close $o
	catch {file rename -force $compar_file $compar_file.old}
	file rename $compar_file.temp $compar_file	
}

proc cg_select_help {} {
set help [file_read $::appdir/lib/cg_multireg.help]
puts [string_change $help [list @BASE@ [get ::base {[info source]}]]]
}

proc cg_multireg {args} {
	if {([llength $args] < 1)} {
		cg_select_help
		exit 1
	}
	foreach {compar_file} $args break
	set files [lrange $args 1 end]
	foreach file $files {
		putslog "Adding $file"
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
