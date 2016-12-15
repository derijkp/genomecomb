#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_multicompar {args} {
	set args [job_init {*}$args]
	set reannot 0
	set regonly 0
	set split 0
	set targetsfile {}
	set targetsfield {}
	cg_options multicompar args {
		-reannot {
			if {$value ni {0 1}} {set value 1; incr pos -1}
			set reannot $value
		}
		-reannotregonly {
			putslog "Also reannot"
			if {$value ni {0 1}} {set value 1; incr pos -1}
			set reannot $value
			set regonly $value
		}
		-split {
			set split [true $value]
		}
		-targetsfile {
			set targetsfile $value
			set targetsfield [lindex [split [file root [file tail $targetsfile]] -] end]
		}
	} compar_file
	set dirs $args
	set compar_file_root [gzroot $compar_file]
	unset -nocomplain a
	if {[file exists $compar_file]} {
		foreach name [cg select -n $compar_file] {
			set a($name) $compar_file
		}
	}
	set files {}
	foreach dir $dirs {
		set dir [file_absolute $dir]
		if {[file isdir $dir]} {
			# a directory is given, look for the variant file(s)
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
			set samples [cg select -n $dir]
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
	# merge variants
	# todo: check for concurrency
	set workdir [indexdir_filewrite $compar_file multicompar]
	# should take into account existing instead of deleting and starting all over -> not now
	if {[file exists $workdir]} {file delete -force $workdir}
	file mkdir $workdir
	job_logdir $workdir/log_jobs
	set multi_merge_num 0
	if {[file exists $compar_file]} {
		file rename -force $compar_file $compar_file.old
		set allfiles [list_concat $compar_file.old $files]
	} else {
		set allfiles $files
	}
	# for calculating the varlines needed, we can treat targetsfile as just another variant file
	set files $allfiles
	if {$targetsfile ne ""} {lappend files $targetsfile}
	multi_merge_job $workdir/vars.tsv $files -split $split
	#
	# make todofile used by multi_join to join all separate files together
	set todofile $workdir/multitodo.txt
	set ftodo [open $todofile w]
	# line with: filename chrpos beginpos endpos typepos refpos altpos todomax todoseqpos ?numberof data positions?
	puts $ftodo [join [list $workdir/vars.tsv 0 1 2 3 4 5 5 0] \t]
	# line with: data positions that must go in the multicompar
	puts $ftodo ""
	set reannotheader {chromosome begin end type ref alt}
	foreach file $allfiles {
		set f [gzopen $file]
		set header [tsv_open $f]
		gzclose $f
		set file [file_absolute $file]
		set dir [file dir $file]
		set filebase [file root [file tail [gzroot $file]]]
		set basicposs [tsv_basicfields $header 6 $file]
		set samples [samples $header]
		if {[llength $samples] > 0} {
			set keepposs [list_find -glob $header *-*]
			# in a multicompar file, sequenced columns should already be present
			# put a value not -1 or -2, so no new one will be added
			set seqpos 0
			lappend reannotheader {*}[list_sub $header $keepposs]
		} else {
			set sample [sourcename $filebase]
			set seqpos [lsearch $header sequenced]
			set keepfields [list_sub $header -exclude $basicposs]
			if {[string match fannotvar-* [file tail $file]]} {
				set mergefields {xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex hasCodingRegion impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence pfam}
				set keepfields [list_lremove $keepfields $mergefields]
			}
			set keepfields [list_remove $keepfields $targetsfield]
			set keepposs [list_cor $header $keepfields]
			if {$seqpos == -1} {
				lappend reannotheader sequenced-$sample
			}
			foreach field $keepfields {
				lappend reannotheader ${field}-$sample
			}
		}
		set max [lmath_max [list_concat $keepposs $basicposs]]
		# line with: filename chrpos beginpos endpos typepos refpos altpos maxcoltoread todoseqpos sizeofnextline
		puts $ftodo [join [list $file {*}$basicposs $max $seqpos [llength $keepposs]] \t]
		# line with: data positioons that must go in the multicompar
		puts $ftodo [join $keepposs \t]
	}
	set len [expr {[llength $allfiles]+1}]
	if {$targetsfile ne ""} {
		set f [open $targetsfile]
		set header [tsv_open $f]
		close $f
		lappend reannotheader $targetsfield
		set file [file_absolute $targetsfile]
		set basicposs [tsv_basicfields $header 6 $file]
		# if seqpos != -1, no new calculated sequence column will be added
		# -2 is used to indicate to put an empty field instead of a ? for no match
		set seqpos -2
		set keepposs [lsearch $header name]
		set max [lmath_max [list_concat $keepposs $basicposs]]
		puts $ftodo [join [list $file {*}$basicposs $max $seqpos [llength $keepposs]] \t]
		puts $ftodo [join $keepposs \t]
		incr len
	}
	close $ftodo
	job multi_join -deps [list $todofile] -vars {split len reannotheader} -targets $compar_file_root -code {
		set tempfile [filetemp $target]
		file_write $tempfile [join $reannotheader \t]\n
		exec multi_join $dep1 $len $split >> $tempfile
		file rename -force $tempfile $target
	}
	if {$reannot} {
		job multicompar_reannot -deps {$compar_file_root} -vars {regonly} -targets {$compar_file_root $compar_file_root.reannot} -code {
			putslog "Reannotating $dep"
			multicompar_reannot $dep 0 $regonly
			file_write $dep.reannot [timestamp]
		}
	}
	job_wait
}
