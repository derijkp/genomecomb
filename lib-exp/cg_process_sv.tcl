#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc process_sv {cgdir dir dbdir {force 0}} {
	job_logdir $dir/log_jobs
	set keepdir [pwd]
	set cgdir [file normalize $cgdir]
	set dir [file normalize $dir]
	set dbdir [file normalize $dbdir]
	set name [file tail $dir]
	file mkdir $dir/sv.temp
	cd $dir/sv.temp

	if {$force || ![file exists $dir/sv/${name}_map2sv_sort_FINISHED]} {
		cg map2sv $cgdir $dir/sv.temp/$name >@ stdout
		file rename {*}[glob $dir/sv.temp/*] $dir/sv
		file delete $dir/sv.temp
	}
	set resultfiles {}
	cd $dir/sv
	set files [gzfiles $name-*-paired.tsv]
	foreach file [lsort -dict $files] {
		set root [file root [gzroot $file]]
		job svindex-$root -deps $file -targets $root.tsv.end1_index -code {
			puts "Indexing $dep"
			cg tsv_index end1 $dep
		}
		job svinfo-$root -deps $file -targets $root.tsv.numinfo -code {
			puts "Info on $dep"
			cg svinfo $dep
		}
		job svfind-$root -deps {$file $root.tsv.numinfo} -targets $root-sv.tsv -vars {dbdir} -code {
			puts "svfind $dep"
			cg svfind $dep [lindex [glob $dbdir/reg_*_simpleRepeat.tsv] 0]
		}
		lappend resultfiles $root-sv.tsv
#		if {$force || [file extension $file] ne ".rz"} {
#			puts "razipping $file"
#			cg_razip $file
#		}
	}
	job svfind-$root -deps $resultfiles -targets $dir/sv-$name.tsv -vars {dbdir} -code {
		cg cat {*}$deps > $target.temp
		file rename $target.temp $target
		putslog "Done: finished finding sv in $dir"
	}
	cd $keepdir
}

proc cg_process_sv {args} {
	global scriptname action
	set args [job_args $args]
	if {([llength $args] < 2) || ([llength $args] > 3)} {
		puts stderr "format is: $scriptname $action sampledir destdir dbdir ?force?"
		puts stderr " - processes sv for one sample directory."
		puts stderr " - By default, only files that are not present already will be created."
		puts stderr " -When force is given as a parameter, everything will be recalculated and overwritten."
		exit 1
	}
	foreach {dir destdir dbdir force} $args break
	switch $force {
		force {set force 1}
		"" {set force 0}
		default {error "unrecognized option $force"}
	}
	process_sv $dir $destdir $dbdir $force
	job_wait
}
