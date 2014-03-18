#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multi_merge_job {varsfile files {split 1}} {
	global multi_merge_num
	upvar job_logdir job_logdir
	set newfiles {}
	if {[llength $files] < 2} {
		mklink [lindex $files 0] $varsfile
		return
	} 
	foreach {file1 file2} $files {
		if {$file2 eq ""} continue
		incr multi_merge_num
		job multi_merge-$multi_merge_num -deps [list $file1 $file2] -vars split -targets $varsfile.$multi_merge_num -code {
			set f [gzopen $dep1]
			set header [tsv_open $f]
			set poss1 [tsv_basicfields $header 6 1]
			close $f
			set f [gzopen $dep2]
			set header [tsv_open $f]
			set poss2 [tsv_basicfields $header 6 1]
			close $f
			exec multi_merge $dep1 {*}$poss1 $dep2 {*}$poss2 $split > $target.temp
			file rename -force $target.temp $target
		}
		lappend newfiles $varsfile.$multi_merge_num
	}
	if {$file2 eq ""} {
		lappend newfiles $file1
	}
	if {[llength $newfiles] > 1} {
		multi_merge_job $varsfile $newfiles $split
	} else {
		incr multi_merge_num
		job multi_merge-$multi_merge_num -deps $newfiles -vars split -targets $varsfile -code {
			file rename -force $dep1 $target
		}
	}
	
}

proc cg_multicompar {args} {
	set args [job_init -silent {*}$args]
	set reannot 0
	set regonly 0
	set split 0
	set listfields {}
	set pos 0
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-reannot {
				putslog "Also reannot"
				set reannot 1
			}
			-reannotregonly {
				putslog "Also reannot"
				set reannot 1
				set regonly 1
			}
			-split {
				incr pos
				set split [true [lindex $args $pos]]
			}
			-listfields {
				incr pos
				set listfields [lindex $args $pos]
			}
			-- break
			default {
				break
			}
		}
		incr pos
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 1)} {
		puts "Wrong number of arguments"
		errorformat multicompar
		exit 1
	}
	foreach {compar_file} $args break
	unset -nocomplain a
	if {[file exists $compar_file]} {
		foreach name [cg select -n $compar_file] {
			set a($name) $compar_file
		}
	}
	file mkdir $compar_file.index
	set workdir $compar_file.index/multicompar
	# should take into account existing instead of deleting and starting all over -> not now
	if {[file exists $workdir]} {file delete -force $workdir}
	file mkdir $workdir
	set dirs [lrange $args 1 end]
	set files {}
	foreach dir $dirs {
		set dir [file_absolute $dir]
		if {[file isdir $dir]} {
			set name [file tail $dir]
			if {[file exists $dir/fannotvar-$name.tsv]} {
				set file [gzfile $dir/fannotvar-$name.tsv]
				lappend files $file
				if {[info exists a($name)]} {
					error "sample \"$name\" in \"$file\" was already in \"$a($name)\""
				}
				set a($name) $file
			} else {
				foreach file [gzfiles $dir/var-*-$name.tsv] {
					set base [file root [file tail [gzroot $file]]]
					set file [gzfile $file]
					set name [sourcename $base]
					lappend files $file
					if {[info exists a($name)]} {
						error "sample \"$name\" in \"$file\" was already in \"$a($name)\""
					}
					set a($name) $file
				}
			}
		} elseif {[llength [cg select -n $dir]]} {
			set header [cg select -h $dir]
			set samples [samples $header]
			set file $dir
			lappend files $file
			foreach name $samples {
				if {[info exists a($name)]} {
					error "sample \"$name\" in \"$file\" was already in \"$a($name)\""
				}
				set a($name) $file
			}
		} else {
			set base [file root [file tail [gzroot $dir]]]
			set name [sourcename $base]
			set file [gzfile $dir]
			lappend files $file
			if {[info exists a($name)]} {
				error "sample \"$name\" in \"$file\" was already in \"$a($name)\""
			}
			set a($name) $file
		}
	}
	job_logdir $workdir/log_jobs
	# merge variants
	# todo: check for concurrency
	set multi_merge_num 0
	if {[file exists $compar_file]} {
		file rename -force $compar_file $compar_file.old
		set allfiles [list_concat $compar_file.old $files]
	} else {
		set allfiles $files
	}
	multi_merge_job $workdir/vars.tsv $allfiles $split
	set todofile $workdir/multitodo.txt
	set ftodo [open $todofile w]
	puts $ftodo [join [list $workdir/vars.tsv 0 1 2 3 4 5 5 0] \t]
	puts $ftodo ""
	set mergefields {xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex hasCodingRegion impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence pfam}
	set reannotheader {chromosome begin end type ref alt}
	foreach file $allfiles {
		set f [open $file]
		set header [tsv_open $f]
		close $f
		set file [file_absolute $file]
		set dir [file dir $file]
		set filebase [file root [file tail [gzroot $file]]]
		set basicposs [tsv_basicfields $header]
		set samples [samples $header]
		if {[llength $samples] > 0} {
			set keepposs [list_find -glob $header *-*]
			set seqpos -2
			lappend reannotheader {*}[list_sub $header $keepposs]
		} else {
			set sample [sourcename $filebase]
			set seqpos [lsearch $header sequenced]
			set keepfields [list_sub $header -exclude $basicposs]
			if {[string match fannotvar-* [file tail $file]]} {
				set keepfields [list_lremove $keepfields $mergefields]
			}
			set keepposs [list_cor $header $keepfields]
			if {$seqpos == -1} {
				lappend reannotheader sequenced-$sample
			}
			foreach field $keepfields {
				lappend reannotheader ${field}-$sample
			}
		}
		set max [lmath_max [list_concat $keepposs $basicposs]]
		puts $ftodo [join [list $file {*}$basicposs $max $seqpos [llength $keepposs]] \t]
		puts $ftodo [join $keepposs \t]
	}
	close $ftodo
	set len [expr {[llength $allfiles]+1}]
	job multi_join -deps [list $todofile {*}$allfiles] -vars {split len reannotheader} -targets $compar_file -code {
		file_write $target.temp [join $reannotheader \t]\n
		exec multi_join $dep1 $len $split >> $target.temp
		file rename -force $target.temp $target
	}
	if {$reannot} {
		job multicompar_reannot -deps {$compar_file} -vars {regonly} -targets {$compar_file $compar_file.reannot} -code {
			putslog "Reannotating $dep"
			multicompar_reannot $dep 0 $regonly
			file_write $dep.reannot [timestamp]
		}
	}
	job_wait
}

