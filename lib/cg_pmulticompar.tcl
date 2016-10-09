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

proc pmulticompar_job {compar_file dirs {regonly 0} {split 1} {targetsfile {}} {erroronduplicates 0} {skipincomplete 0}} {
# putsvars compar_file dirs regonly split targetsfile erroronduplicates
	if {[jobfileexists $compar_file]} {
		set dirs [list $compar_file {*}$dirs]
	}
	set compar_file [file_absolute $compar_file]
	set basedir [file dir $compar_file]
	set compar_file_root [gzroot $compar_file]
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
			if {[file exists $dir]} {
				set filesamples [cg select -n $dir]
			} else {
				set filesamples {}
			}
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
	set workdir $compar_file.index/multicompar
	# should take into account existing instead of deleting and starting all over -> not now
	if {[file exists $workdir]} {file delete -force $workdir}
	file delete -force $workdir
	file mkdir $workdir
	job_logdir $workdir/log_jobs
	set multi_merge_num 0
	if {[file exists $compar_file]} {
		# file rename -force $compar_file $compar_file.old
		set allfiles [list_concat $compar_file $files]
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
		set varallfile $sampledir/varall-$sample.tsv
		set target $workdir/avars-$sample.tsv
		lappend pastefiles $target
		if {![jobfileexists $sregfile] && ![jobfileexists $varallfile]} {
			set msg "no sorted region file ($sregfile) or varallfile ($varallfile) found: not properly processed sample"
			if {!$skipincomplete} {
				error $msg
			} else {
				putslog "warning: $msg"
			}
		}
		job multicompar_addvars-$sample -deps {$allvarsfile $samplevarsfile ($sregfile) ($varallfile)} -targets {$target} \
		  -vars {allvarsfile samplevarsfile sregfile varallfile sample split} -code {
			set vf [gzopen $allvarsfile]
			set vheader [tsv_open $vf]
			close $vf
			if {$vheader ne "chromosome begin end type ref alt"} {error "internal error: index vars.tsv in wrong format"}
			set f [gzopen $samplevarsfile]
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
			set sregfile [lindex [gzfiles $sregfile] 0]
			set varallfile [lindex [gzfiles $varallfile] 0]
			# puts [list ../bin/multicompar_addvars_simple $allvarsfile $samplevarsfile $split $sregfile $varallfile {*}$keepposs]
			exec multicompar_addvars_simple $allvarsfile $samplevarsfile $split $sregfile $varallfile {*}$keepposs >> $target.temp
			file rename $target.temp $target
		}
	}
	if {$targetsfile ne ""} {
		# add targetsfile annotation
		set target $workdir/targets_annot.tsv
		lappend pastefiles $target
		job multicompar_targets -deps {$allvarsfile $targetsfile} -targets {$target} \
		  -vars {allvarsfile targetsfile split} -code {
			set targetsfield [lindex [split [file root [file tail $targetsfile]] -] end]
			set f [gzopen $targetsfile]
			set header [tsv_open $f]
			gzclose $f
			set tempdbposs [tsv_basicfields $header 6 0]
			set dbposs [lrange $tempdbposs 0 2]
			set type2pos [lindex $tempdbposs 3]
			if {$type2pos == -1} {
				error "$targetsfile has no type field"
			}
			set alt2pos [lindex $tempdbposs 5]
			if {$alt2pos == -1} {
				error "$targetsfile has no alt field"
			}
			set keeppos [lsearch $header name]
			if {$keeppos == -1} {set keeppos {}}
			file_write $target.temp $targetsfield\n
			# puts [list ../bin/var_annot $allvarsfile 0 1 2 3 5 $targetsfile {*}$dbposs $type2pos $alt2pos "" {*}$keeppos]
			exec var_annot $allvarsfile 0 1 2 3 5 $targetsfile {*}$dbposs $type2pos $alt2pos "" {*}$keeppos >> $target.temp 2>@ stderr
			file rename $target.temp $target
		}
	}
	#
	# paste all together for final multicompar
	tsv_paste_job $compar_file $pastefiles -forcepaste 1 -endcommand [list file delete {*}$pastefiles]
}

proc cg_pmulticompar {args} {
	set args [job_init -silent {*}$args]
	set regonly 0
	set split 0
	set erroronduplicates 0
	set targetsfile {}
	set skipincomplete 1
	cg_options pmulticompar args {
		-r - -reannotregonly - --reannotregonly {
			putslog "Reannot reg only"
			set regonly $value
		}
		-s - -split - --split {
			set split $value
		}
		-e - --erroronduplicates {
			set erroronduplicates $value
		}
		-t - -targetsfile - --targetsfile {
			set targetsfile $value
		}
		-i - -skipincomplete - --skipincomplete {
			set skipincomplete $value
		}
	} 1
	foreach {compar_file} $args break
	set dirs [lrange $args 1 end]
	pmulticompar_job $compar_file $dirs $regonly $split $targetsfile $erroronduplicates $skipincomplete
	job_wait
}

