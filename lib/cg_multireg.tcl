#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multireg_job {compar_file regfiles {limitreg {}}} {
	set compar_file [file_absolute $compar_file]
	job_logdir [gzroot $compar_file].index/log_jobs
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
	if {$limitreg ne ""} {lappend files $limitreg}
	set len [llength $files]
	if {$len <= $maxfiles} {
		set target $compar_file
		job multireg-[file tail $compar_file] -force $jobforce -deps $files -targets {$target} -vars {isreg limitreg} -code {
			if {$limitreg ne ""} {
				set templist {}
				foreach file $deps {
					set tempfile [tempdir]/[file tail $file]
					lappend templist $tempfile
					exec cg regselect -o $tempfile $file $limitreg
				}
				set deps $templist
			}
			set compress [compresspipe $target]
			set todo [list_merge $deps $isreg]
			set temp [filetemp_ext $target]
			# puts [list ../bin/multireg {*}$todo]
			exec multireg {*}$todo {*}$compress > $temp 2>@ stderr
			if {$limitreg ne ""} {
				set temp2 [filetemp_ext $target]
				set l [file root [file tail [gzroot $limitreg]]]
				cg select -overwrite 1 -f [list_remove [cg select -h $temp] $l] -q "\$$l == 1" $temp $temp2
				file rename -force $temp2 $target
				file delete $temp
			} else {
				file rename -force $temp $target
			}
			if {$compress ne ""} {cg_zindex $target}
		}
		return
	}
	set workdir [gzroot $compar_file].index/multicompar
	file mkdir $workdir
	catch {file delete {*}[glob -nocomplain $workdir/multireg.temp*]}
	set todo $files
	set todoisreg $isreg
	set delete 0
	set num 1
	while 1 {
		if {$len <= $maxfiles} {
			set target $compar_file
			job multireg-[file tail $compar_file] -force $jobforce -deps $todo -targets {
				$target
			} -vars {
				todoisreg delete workdir limitreg
			} -code {
				set compress [compresspipe $target]
				set todo [list_merge $deps $todoisreg]
				# puts [list ../bin/multireg {*}$todo]
				set temp [filetemp_ext $target]
				exec multireg {*}$todo {*}$compress > $temp 2>@ stderr
				if {$limitreg ne ""} {
					set temp2 [filetemp_ext $target]
					set l [file root [file tail [gzroot $limitreg]]]
					cg select -overwrite 1 -f [list_remove [cg select -h $temp] $l] -q "\$$l == 1" $temp $temp2
					file rename -force $temp2 $target
				} else {
					file rename -force $temp $target
				}
				if {$compress ne ""} {cg_zindex $target}
				if {$delete} {file delete {*}$deps}
				file delete $workdir
			}
			break
		}
		set pos 0
		set newtodo {}
		while {$pos < $len} {
			# go over sets (size maxfiles) and merge to temporary files that will be merged in the final merge
			set deps [lrange $todo $pos [expr {$pos+$maxfiles-1}]]
			set partisreg [lrange $todoisreg $pos [expr {$pos+$maxfiles-1}]]
			if {[llength $deps] > 1} {
				set target $workdir/multireg.temp$num.zst
				incr num
				lappend newtodo $target
				lappend newisreg 1
				job multireg-[file tail $target] -deps $deps -targets {
					$target
				} -vars {
					partisreg delete limitreg
				} -code {
					# puts [list ../bin/multireg {*}$deps]
					if {!$delete && $limitreg ne ""} {
						set templist {}
						set tempdir [tempdir]
						foreach file $deps {
							set tempfile $tempdir/[file tail $file]
							lappend templist $tempfile
							exec cg regselect -o $tempfile $file $limitreg
						}
						set deps $templist
					}
					set todo [list_merge $deps $partisreg]
					# puts [list ../bin/multireg {*}$todo]
					set temp [filetemp $target]
					exec multireg {*}$todo {*}[compresspipe $target -1] > $temp 2>@ stderr
					file rename -force $temp $target
					if {$delete} {file delete {*}$deps}
				}
			} else {
				set dep [lindex $deps 0]
				set target $workdir/[file tail $dep]
				lappend newtodo $target
				lappend newisreg [lindex $partisreg 0]
				job multireg-[file tail $target] -deps $deps -targets {
					$target
				} -vars {
					partisreg delete limitreg
				} -code {
					if {$limitreg ne ""} {
						exec cg regselect -o $target.temp $dep $limitreg
					} elseif {$delete} {
						file rename $dep $target.temp
					} else {
						mklink $dep $target.temp
					}
					file rename -force $target.temp $target
				}
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
	set limitreg {}
	cg_options multireg args {
		-m - -maxopenfiles {
			set ::maxopenfiles [expr {$value - 4}]
		}
		-limitreg {
			set limitreg $value
		}
	} compar_file 2
	multireg_job $compar_file $args $limitreg
}
