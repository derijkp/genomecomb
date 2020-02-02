proc process_multicompar_job {args} {
	set keepargs $args
	set dbdir {}
	set dbfiles {}
	set realign 1
	set cleanup 1
	set paired 1
	set conv_nextseq 0
	set skipincomplete 1
	set targetvarsfile {}
	set addtargets 0
	set threads 1
	set distrreg 0
	set keepfields *
	set limitreg {}
	cg_options process_multicompar args {
		-dbdir {
			set dbdir $value
		}
		-split {
			set split $value
		}
		-dbfile {
			set file [gzfile $value]
			if {![file exists $file]} {error "dbfile $value does not exists"}
			lappend dbfiles [file_absolute $file]
		}
		-dbfiles {
			foreach v $value {
				set file [gzfile $v]
				if {![file exists $file]} {error "dbfile $v does not exist"}
				lappend dbfiles [file_absolute $file]
			}
		}
		-skipincomplete {
			set skipincomplete $value
		}
		-targetvarsfile {
			set addtargets 1
			if {$targetvarsfile ne "" && ![file exists $value]} {error "targetvarsfile $value does not exists"}
			set targetvarsfile $value
		}
		-varfiles {
			set varfiles {}
			foreach file $value {
				lappend varfiles [file_absolute $file]
			}
		}
		-svfiles {
			set svfiles $value
		}
		-experiment {
			set experiment $value
		}
		-threads {
			set threads $value
			# not used yet
		}
		-distrreg {
			set distrreg [distrreg_checkvalue $value]
		}
		-keepfields {
			set keepfields $value
		}
		-limitreg {
			set limitreg $value
		}
		-c - -cleanup {
			set cleanup $value
		}
		-m - -maxopenfiles {
			set ::maxopenfiles [expr {$value - 4}]
		}
	} {destdir dbdir} 1
	set dbfiles [list_remove [list_remdup $dbfiles] {}]
	set destdir [file_absolute $destdir]
	file mkdir $destdir/compar
	if {![info exists experiment]} {
		set experiment [file tail $destdir]
	}
	# check projectinfo
	projectinfo $destdir dbdir {split 1}
	set dbdir [dbdir $dbdir]
	set refseq [glob $dbdir/genome_*.ifas]
	# logfile
	# -------
	set cmdline [list cg process_multicompar]
	foreach option {
		dbdir split dbfile dbfiles skipincomplete targetvarsfile varfiles experiment cleanup maxopenfiles
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline $destdir
	upvar job_logdir job_logdir
	job_logfile $destdir/process_multicompar_[file tail $destdir] $destdir $cmdline \
		{*}[versions dbdir gnusort8 zst os]
	# analysis info
	# -------------
	info_analysis_file $destdir/compar/info_analysis.tsv {} \
		{dbdir split dbfiles targetvarsfile ::maxopenfiles} \
		{genomecomb dbdir gnusort8 zst os} \
		command [list cg process_multicompar {*}$keepargs]

	# get files todo
	# ----------------
	if {[file exists $destdir/samples]} {
		set sampledir $destdir/samples
	} else {
		set sampledir $destdir
	}
	if {![info exists varfiles]} {
		set varfiles [bsort [jobglob ${sampledir}/*/var-*.tsv]]
	} else {
		set tempvarfiles {}
		foreach varfile $varfiles {
			set temp [jobglob $sampledir/*/$varfile]
			if {$temp eq ""} {
				set temp [jobglob $varfile]
			}
			if {$temp eq ""} {
				set temp [jobglob $sampledir/*/[file tail $varfile]]
			}
			if {$temp eq ""} {error "varfile $varfile not found"}
			lappend tempvarfiles {*}[bsort $temp]
		}
		set varfiles $tempvarfiles
	}
	if {![info exists svfiles]} {
		set svfiles [bsort [jobglob ${sampledir}/*/sv-*.tsv]]
	} else {
		set tempsvfiles {}
		foreach svfile $svfiles {
			if {![jobfileexists $svfile]} {
				set temp [jobglob $sampledir/*/$svfile]
				if {$temp eq ""} {
					set temp [jobglob $svfile]
				}
				if {$temp eq ""} {error "svfile $svfile not found"}
				lappend tempsvfiles {*}[bsort $temp]
			} else {
				lappend tempsvfiles $svfile
			}
		}
		set svfiles $tempsvfiles
	}
	# set up job
	job_logdir $destdir/log_jobs
	set keeppwd [pwd]
	# change to workdir
	set keeppwd [pwd]
	cd $destdir
	file mkdir compar
	#
	# multicompar
	# -----------
	# todo: no check if targetsfile was actually used
	if {![llength $varfiles]} {
		putslog "No qualifying varfiles: not making var compar"
	} else {
		putslog "Making/updating multicompar in $destdir/compar/compar-$experiment.tsv"
		putslog "Finding samples"
		set compar_file [gzfile compar/compar-$experiment.tsv]
		if {[catch {cg select -a $compar_file} done]} {set done {}} else {set done [split $done \n]}
		if {[file exists $compar_file]} {set mtime [file mtime $compar_file]} else {set mtime 0}
		set stilltodo {}
		foreach varfile $varfiles {
			if {![file exists $varfile] || [file mtime $varfile] > $mtime} {
				putslog "redo all: $varfile is newer than $compar_file"
				set stilltodo $varfiles
				if {[file exists $compar_file]} {file rename -force -- $compar_file $compar_file.old}
				break
			}
			set analysis [file_analysis $varfile]
			if {![inlist $done $analysis]} {
				putslog "Still todo: $analysis"
				lappend stilltodo $varfile
			}
		}
		putslog "Samples: [llength $varfiles] todo, [llength $done] already done, [llength $stilltodo] to add"
		if {[llength $stilltodo]} {
			putslog "Samples to add: $stilltodo"
			putslog "Starting multicompar"
			set compar_file [gzroot $compar_file].zst
			# pmulticompar_job $compar_file $stilltodo 0 $split $targetvarsfile 0 $skipincomplete
			pmulticompar_job -reannotregonly 0 -split $split -limitreg $limitreg -targetvarsfile $targetvarsfile -erroronduplicates 0 -skipincomplete $skipincomplete -keepfields $keepfields $compar_file {*}$stilltodo
		} else {
			putslog "skipping multicompar (no update needed)"
		}
		# annotate multicompar
		# --------------------
		putslog "Starting annotation"
		cg_annotate_job -distrreg $distrreg $compar_file compar/annot_compar-$experiment.tsv.zst $dbdir {*}$dbfiles
		job indexannotcompar-$experiment \
		-deps {compar/annot_compar-$experiment.tsv} \
		-targets {compar/annot_compar-$experiment.tsv.index/info.tsv} -vars dbdir -code {
			cg index -colinfo $dep
		}
		#
		# multi sreg
		# ----------
		putslog "Starting multisreg"
		set regfiles {}
		foreach varfile $varfiles {
			set analysis [file_analysis $varfile]
			set name [lindex [split $analysis -] end]
			set file $sampledir/$name/sreg-$analysis.tsv
			if {![jobfileexists $file]} {
				if {!$skipincomplete} {
					error "file $file not found"
				} else {
					putslog "warning: file $file not found"
					continue
				}
			}
			lappend regfiles $file
		}
		multireg_job compar/sreg-$experiment.tsv.zst $regfiles $limitreg

	}
	#
	# sv
	# ----
	if {![llength $svfiles]} {
		putslog "No qualifying svfiles: not making sv compar"
	} else {
		putslog "Starting sv"
		set svcompar_file compar/sv-${experiment}.tsv.zst
		set target $svcompar_file
		testmultitarget $target $svfiles
		job sv_multicompar -optional 1 -deps $svfiles -targets {$target} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -a $target.temp]
			} else {
				set done {}
			}
			set todo {}
			foreach file $deps {
				regexp {sv-(.*)\.tsv*} [file tail $file] temp name
				if {![inlist $done $name]} {
					lappend todo $file
				}
			}
			if {[llength $done]} {
				puts "Multisv already present: $done"
			}
			if {[llength $todo]} {
				cg svmulticompar -overlap 40 $target.temp {*}$todo
			}
			cg_zst $target.temp
			file rename -force -- $target.temp.zst $target
		}
		# annotate svmulticompar
		# --------------------
		putslog "Starting annotation"
		cg_annotate_job -distrreg $distrreg $svcompar_file compar/annot_sv-$experiment.tsv.zst $dbdir {*}$dbfiles
		job indexannotcompar-$experiment -deps {
			compar/annot_sv-$experiment.tsv
		} -targets {
			compar/annot_sv-$experiment.tsv.index/info.tsv
		} -vars dbdir -code {
			cg index -colinfo $dep
		}
	}
	#
	# cgsv
	# ----
	putslog "Starting cgsv"
	set files [jobglob $sampledir/*/cgsv-*.tsv]
	if {[llength $files]} {
		set target compar/cgsv-${experiment}.tsv.zst
		testmultitarget $target $files
		job cgsv_multicompar -optional 1 -deps $files -targets {$target} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -a $target.temp]
			} else {
				set done {}
			}
			set todo {}
			foreach file $deps {
				regexp {cgsv-(.*)\.tsv*} [file tail $file] temp name
				if {![inlist $done $name]} {
					lappend todo $file
				}
			}
			if {[llength $done]} {
				puts "Multicgsv already present: $done"
			}
			if {[llength $todo]} {
				cg svmulticompar $target.temp {*}$todo
			}
			cg_zst $target.temp
			file rename -force -- $target.temp.zst $target
		}
		job cgsv_annotate -optional 1 \
		-deps {compar/cgsv-$experiment.tsv} \
		-targets {compar/annot_cgsv-$experiment.tsv.zst} -vars {refseqdir build} -code {
			cg annotate $dep $target $refseqdir/$build
		}
		job cgsv_annotate_index -optional 1 \
		-deps {compar/annot_cgsv-$experiment.tsv.zst} \
		-targets {compar/annot_cgsv-$experiment.tsv.index/info.tsv} -code {
			cg index $dep
		}
	}
	# cgcnv
	# ----
	putslog "Starting cgcnv"
	set files [jobglob $sampledir/*/cgcnv-*.tsv]
	if {[llength $files]} {
		set target compar/cgcnv-${experiment}.tsv.zst
		set names [list_regsub {.*/cgcnv-(.*)\.tsv.*} $files {\1}]
		testmultitarget $target $files
		job cgcnv_multicompar -optional 1 -deps $files -targets {compar/cgcnv-${experiment}.tsv.zst} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -a $target.temp]
			} else {
				set done {}
			}
			set todo {}
			foreach file $deps {
				regexp {cgcnv-(.*)\.tsv*} [file tail $file] temp name
				if {![inlist $done $name]} {
					lappend todo $file
				}
			}
			if {[llength $done]} {
				puts "Multicgcnv already present: $done"
			}
			if {[llength $todo]} {
				cg svmulticompar $target.temp {*}$todo
			}
			cg_zst $target.temp
			file rename -force -- $target.temp.zst $target
		}
		job cgcnv_annotate -optional 1 \
		-deps {compar/cgcnv-$experiment.tsv.zst} \
		-targets {compar/annot_cgcnv-$experiment.tsv.zst} -vars {refseqdir build} -code {
			cg annotate $dep $target $refseqdir/$build
		}
		job cgcnv_annotate_index -optional 1 \
		-deps {compar/annot_cgcnv-$experiment.tsv.zst} \
		-targets {compar/annot_cgcnv-$experiment.tsv.index/info.tsv} -code {
			cg index -colinfo $dep
		}
	}
	cd $keeppwd
}

proc cg_process_multicompar {args} {
	putslog "Running [list cg process_multicompar {*}$args]"
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_multicompar
	}
	process_multicompar_job {*}$args
	job_wait
}

