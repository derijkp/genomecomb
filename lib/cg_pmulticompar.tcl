#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc pmulticompar_findsample {basedir sample args} {
	if {![llength $args]} {set args [list {}]}
	set sampledir [lindex [split $sample -] end]
	foreach startdir [list $basedir/samples $basedir [file dir $basedir]/samples [file dir $basedir]] {
		foreach pattern $args {
			foreach usesample [list $sample $sampledir {}] {
				set test [file join $startdir $usesample $pattern]
				if {[jobfileexists $test]} {
					return $test
				}
			}
		}
	}
	return {}
}

proc pmulticompar_job {newcompar_file dirs {regonly 0} {split 1} {targetsfile {}} {erroronduplicates 0}} {

	set newcompar_file [file_absolute $newcompar_file]
	set basedir [file dir $newcompar_file]
	set targetsfield [lindex [split [file root [file tail $targetsfile]] -] end]
	set newcompar_file_root [gzroot $newcompar_file]
	unset -nocomplain samplesa
	set files {}
	set samples {}
	foreach dir $dirs {
		set dir [file_absolute $dir]
		if {![jobfileexists $dir]} {error "$dir does not exist"}
		if {[file isdir $dir]} {
			# a directory is given, look for the variant file(s)
			set sample [file tail $dir]
			set file [lindex [jobglob $dir/fannotvar-$sample.tsv] 0]
			if {$file ne ""} {
				lappend files $file
				if {[info exists samplesa($sample)]} {
					puts stderr "sample \"$sample\" in \"$file\" was already in \"$samplesa($sample)\""
					if {$erroronduplicates} {exit 1} else continue
				}
				set samplesa($sample) $file
				lappend samples $sample
			} else {
				set files [jobglob $dir/var-*-$sample.tsv]
				if {![llength $files]} {error "no variant files found in dir $dir"}
				foreach file $files {
					set base [file root [file tail [gzroot $file]]]
					set file [gzfile $file]
					set sample [sourcename $base]
					if {[info exists samplesa($sample)]} {
						puts stderr "sample \"$sample\" in \"$file\" was already in \"$samplesa($sample)\""
						if {$erroronduplicates} {exit 1} else continue
					}
					lappend files $file
					set samplesa($sample) $file
					lappend samples $sample
				}
			}
		} else {
			set filesamples [cg select -n $dir]
			if {[llength $filesamples]} {
				# is a multicompar file
				set file $dir
				set newsample 0
				foreach sample $filesamples {
					if {[info exists samplesa($sample)]} {
						puts stderr "sample \"$sample\" in \"$file\" was already in \"$samplesa($sample)\""
						if {$erroronduplicates} {exit 1} else continue
					}
					# find ori variant file (and thus sampledir)
					set samplefile [pmulticompar_findsample $basedir $sample fannotvar-$sample.tsv var-$sample.tsv var-*-$sample.tsv]
					if {$samplefile eq ""} {
						error "Variant file for sample $sample not found"
					}
					set samplesa($sample) $samplefile
					lappend samples $sample
					set newsample 1
				}
				if {$newsample} {lappend files $file}
			} else {
				set base [file root [file tail [gzroot $dir]]]
				set sample [sourcename $base]
				set file $dir
				if {[info exists samplesa($sample)]} {
					puts stderr "sample \"$sample\" in \"$file\" was already in \"$samplesa($sample)\""
					if {$erroronduplicates} {exit 1} else continue
				}
				lappend files $file
				set samplesa($sample) $file
				lappend samples $sample
			}
		}
	}
	#
	# merge variants
	# todo: check for concurrency
	set workdir [indexdir_filewrite $newcompar_file multicompar]
	# should take into account existing instead of deleting and starting all over -> not now
	if {[file exists $workdir]} {file delete -force $workdir}
	file mkdir $workdir
	job_logdir $workdir/log_jobs
	set multi_merge_num 0
	if {[file exists $newcompar_file]} {
		file rename -force $newcompar_file $newcompar_file.old
		set allfiles [list_concat $newcompar_file.old $files]
	} else {
		set allfiles $files
	}
	# for calculating the varlines needed, we can treat targetsfile as just another variant file
	set files $allfiles
	if {$targetsfile ne ""} {lappend files $targetsfile}
	multi_merge_job $workdir/vars.tsv $files $split

	# 
	# add extra var line to each sample file to get everything in vars.tsv
	set pastefiles [list $workdir/vars.tsv]
	set allvarsfile $workdir/vars.tsv
	foreach sample $samples {
		set samplevarsfile $samplesa($sample)
		set sampledir [file dir $samplesa($sample)]
		set sregfile $sampledir/sreg-$sample.tsv
		set target $workdir/avars-$sample.tsv
		lappend pastefiles $target
		job multicompar_fill-$sample -deps {$allvarsfile $samplevarsfile ($sregfile)} -targets {$target} \
		  -vars {allvarsfile samplevarsfile sregfile sample split} -code {
			set vf [open $allvarsfile]
			set vheader [tsv_open $vf]
			close $vf
			if {$vheader ne "chromosome begin end type ref alt"} {error "internal error: index vars.tsv in wrong format"}
			set f [open $samplevarsfile]
			set header [tsv_open $f]
			close $f
			set basicposs [tsv_basicfields $header 6 0]
			set seqpos [lsearch $header sequenced]
			set zygpos [lsearch $header zyg]
			set a1pos [lsearch $header alleleSeq1]
			set a2pos [lsearch $header alleleSeq2]
			set keepfields [list_sub $header -exclude $basicposs]
			if {[string match fannotvar-* [file tail $samplevarsfile]]} {
				set mergefields {xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex hasCodingRegion impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence pfam}
				set keepfields [list_lremove $keepfields $mergefields]
			}
			set keepfields [list_remove $keepfields sequenced zyg alleleSeq1 alleleSeq2]
			set keepposs [list_cor $header $keepfields]
			if {$seqpos == -1} {
				lappend reannotheader sequenced-$sample
			}
			set newheader [list sequenced-$sample]
			if {$zygpos != -1} {lappend newheader zyg-$sample}
			if {$a1pos != -1} {lappend newheader alleleSeq1-$sample}
			if {$a2pos != -1} {lappend newheader alleleSeq2-$sample}
			foreach field $keepfields {
				lappend newheader ${field}-$sample
			}
			set o [open $target.temp w]
			puts $o [join $newheader \t]
			close $o
			if {![file exists $sregfile]} {set sregfile {}}
			 puts [list ../bin/multicompar_addvars_simple $allvarsfile $samplevarsfile $split $sregfile {*}$keepposs]
			exec multicompar_addvars_simple $allvarsfile $samplevarsfile $split $sregfile {*}$keepposs >> $target.temp
			file rename $target.temp $target
		}
	}
	#
	# paste all together for final multicompar
	tsv_paste_job $newcompar_file $pastefiles
}

proc cg_pmulticompar {args} {
	set args [job_init -silent {*}$args]
	set regonly 0
	set split 0
	set erroronduplicates 0
	set targetsfile {}
	set targetsfield {}
	cg_options pmulticompar args {
		-r - -reannotregonly - --reannotregonly {
			putslog "Reannot reg only"
			set regonly $value
		}
		-s - -split - --split {
			set split $value
		}
		-d - --erroronduplicates {
			set erroronduplicates $value
		}
		-t - -targetsfile - --targetsfile {
			set targetsfile $value
		}
	} 2
	foreach {newcompar_file} $args break
	pmulticompar_job $newcompar_file [lrange $args 1 end] $regonly $split $targetsfile $erroronduplicates
	job_wait
}

