#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc process_sv {cgdir dir dbdir {force 0}} {
	job_logdir $dir/log_jobs
	set keepdir [pwd]
	set cgdir [file_absolute $cgdir]
	set dir [file_absolute $dir]
	set dbdir [file_absolute $dbdir]
	set name [file tail $dir]
	file mkdir $dir/sv.temp
	file mkdir $dir/sv
	cd $dir/sv.temp

	if {$force || ![file exists $dir/sv/${name}_map2sv_sort_FINISHED]} {
		cg map2sv $cgdir $dir/sv.temp/$name >@ stdout
		file rename -force -- {*}[glob $dir/sv.temp/*] $dir/sv
		file delete $dir/sv.temp
	}
	set resultfiles {}
	cd $dir/sv
	set files [gzfiles $name-*-paired.tsv]
	foreach file [bsort $files] {
		set root [file root [gzroot $file]]
		job svindex-$root -deps {$file} -targets {$root.tsv.end1_index} -code {
			puts "Indexing $dep"
			cg tsv_index end1 $dep
		}
		job svinfo-$root -deps {$file} -targets {$root.tsv.numinfo} -code {
			puts "Info on $dep"
			cg svinfo $dep
		}
		job svfind-$root -deps {$file $root.tsv.numinfo} -targets {$root-sv.tsv} -vars {dbdir} -code {
			puts "svfind $dep"
			cg svfind $dep [gzfile $dbdir/reg_*_simpleRepeat.tsv]
		}
		lappend resultfiles $root-sv.tsv
#		if {$force || [file extension $file] ne ".rz"} {
#			puts "razipping $file"
#			cg_razip_job $file
#		}
	}
	job svfind-[file_part $dir end] -deps $resultfiles -targets {$dir/sv/svall-$name.tsv $dir/sv-$name.tsv} -vars {dbdir dir} -code {
		set temptarget [filetemp $target]
		cg cat {*}[bsort $deps] > $temptarget
		file rename -force -- $temptarget $target
		set temptarget [filetemp $target2]
		cg select -overwrite 1 -q {$problems eq "" and $quality > 2} $target $temptarget
		file rename -force -- $temptarget $target2
		putslog "Done: finished finding sv in $dir"
	}
	cd $keepdir
}

proc cg_process_sv {args} {
	global scriptname action
	set args [job_init {*}$args]
	set force 0
	set dbdir {}
	cg_options process_sv args {
		-force {set force $value}
	} {dir destdir dbdir force} 2 4
	switch $force {
		1 - 0 {# keep value}
		force {set force 1}
		"" {set force 0}
		default {error "unrecognized option $force"}
	}
	set dbdir [dbdir $dbdir]
	process_sv $dir $destdir $dbdir $force
	job_wait
}
