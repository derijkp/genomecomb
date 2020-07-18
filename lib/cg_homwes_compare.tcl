#!/bin/sh
# the next line restarts using wish \
exec cg source "$0" ${1+"$@"}

proc findfield {fields pattern} {
	lindex $fields [lsearch -glob $fields $pattern]
}

proc cg_homwes_compare args {
	set mutireg {}
#	set pos 0
#	foreach {key value} $args {
#		switch -- $key {
#			-mutireg {
#				set mutireg $value
#			}
#			default break
#		}
#		incr pos 2
#	}
#	set args [lrange $args $pos end]
	cg_options homwes_compare args {
	} {resultfile} 2
	if {[file exists $resultfile]} {file rename -force -- $resultfile $resultfile.old}
	if {[file exists $resultfile.summary.tsv]} {file rename -force -- $resultfile.summary.tsv $resultfile.summary.tsv.old}
	cg multireg $resultfile {*}$args
	set o [open $resultfile.summary.tsv w]
	puts $o [join {name snps numregions covered common sensitivity specificity fp fn fpvscovered fnvsrefcovered} \t]
	set reffile [lindex $args 0]
	set refname [file tail [file root $reffile]]
	set refcovered [lindex [cg covered $reffile] end]
	unset -nocomplain donea
	foreach file $args {
		set name [file tail [file root $file]]
		if {[info exists donea($name)]} {
			error "duplicate name $name (file $file and $donea($name))"
		}
		putslog "Analysing $name"
		set numregions [lindex [cg select -g all $file] end]
		set covered [lindex [cg covered $file] end]
		set filtered [glob -nocomplain [file root $file].work/$name*-filtered.tsv]
		if {[llength $filtered] == 1} {
			set snps [lindex [cg select -g all [lindex $filtered 0]] end]
		} else {
			set snps NA
		}
		set common [lindex [exec cg select -q "\$$refname == 1 and \$$name == 1" $resultfile | cg covered] end]
		set sensitivity [expr {100.0*$common/$refcovered}]
		set specificity [expr {100.0*$common/$covered}]
		set fp [lindex [exec cg select -q "\$$refname == 0 and \$$name == 1" $resultfile | cg covered] end]
		set fn [lindex [exec cg select -q "\$$refname == 1 and \$$name == 0" $resultfile | cg covered] end]
		set fppercent [expr {100.0*$fp/$covered}]
		set fnpercent [expr {100.0*$fn/$refcovered}]
		puts $o [join [list $name $snps $numregions $covered $common $sensitivity $specificity $fp $fn $fppercent $fnpercent] \t]
		set donea($name) $file
	}
	if {$o ne "stdout"} {close $o}
}
