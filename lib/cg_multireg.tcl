#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multireg_job {compar_file regfiles} {
	set compar_file [file_absolute $compar_file]
	job_logdir $compar_file.index/log_jobs
	set maxfiles [maxopenfiles]
	if {$maxfiles < 2} {set maxfiles 2}
	set fieldsneeded {}
	set files {}
	set todo {}
	if {[file exists $compar_file]} {
		set jobforce 1
		set f [gzopen $compar_file]
		set header [tsv_open $f]
		gzclose $f
		lappend files $compar_file
		lappend isreg 1
	} else {
		set header {}
		set jobforce 0
	}
	foreach file $regfiles {
		set name [file root [file tail [gzroot $file]]]
		if {[inlist $header $name]} {
			putslog "*** Skipping $file: $name already in $compar_file ***"
			continue
		}
		putslog "Adding $file to $compar_file"
		lappend files $file
		lappend isreg 0
		lappend fieldsneeded $name
	}
	if {![llength $files]} return
	set len [llength $files]
	if {$len <= $maxfiles} {
		set target $compar_file
		job multireg-[file tail $compar_file] -force $jobforce -deps $files -targets {$target} -vars {isreg} -code {
			set todo [list_merge $deps $isreg]
			# puts [list ../bin/multireg {*}$todo]
			exec multireg {*}$todo > $target.temp 2>@ stderr
			file rename -force $target.temp $target
		}
		return
	}
	set workdir [indexdir_filewrite $compar_file multireg]
	file mkdir $workdir
	catch {file delete {*}[glob -nocomplain $workdir/multireg.temp*]}
	set todo $files
	set todoisreg $isreg
	set delete 0
	set num 1
	while 1 {
		if {$len <= $maxfiles} {
			set target $compar_file
			job multireg-[file tail $compar_file] -force $jobforce -deps $todo -targets {$target} -vars {todoisreg delete workdir} -code {
				set todo [list_merge $deps $todoisreg]
				# puts [list ../bin/multireg {*}$todo]
				exec multireg {*}$todo > $target.temp 2>@ stderr
				file rename -force $target.temp $target
				# if {$delete} {file delete {*}$deps $workdir}
			}
			break
		}
		set pos 0
		set newtodo {}
		while {$pos < $len} {
			set deps [lrange $todo $pos [expr {$pos+$maxfiles-1}]]
			set partisreg [lrange $todoisreg $pos [expr {$pos+$maxfiles-1}]]
			if {[llength $deps] > 1} {
				set target $workdir/multireg.temp$num
				incr num
				lappend newtodo $target
				lappend newisreg 1
				job multireg-[file tail $target] -deps $deps -targets {$target} -vars {partisreg delete} -code {
					# puts [list ../bin/multireg {*}$deps]
					if {[llength $deps] > 1} {
						set todo [list_merge $deps $partisreg]
						# puts [list ../bin/multireg {*}$todo]
						exec multireg {*}$todo > $target.temp 2>@ stderr
					} elseif {!$delete} {
						mklink $dep $target.temp
					} else {
						file rename $dep $target.temp
					}
					file rename -force $target.temp $target
					#if {$delete} {file delete {*}$deps}
				}
			} else {
				lappend newtodo [lindex $deps 0]
				lappend newisreg [lindex $partisreg 0]
			}
			incr pos $maxfiles
			
		}
		set delete 1
		set todo $newtodo
		set todoisreg $newisreg
		set len [llength $todo]
	}
}

proc cg_multireg {args} {
	set args [job_init {*}$args]
	cg_options multireg args {
		-m - --maxopenfiles {
			set ::maxopenfiles [expr {$value - 4}]
		}
	} compar_file 2
	multireg_job $compar_file $args
}
