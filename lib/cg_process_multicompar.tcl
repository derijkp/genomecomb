proc process_multicompar_job {args} {
	set keepargs $args
	set dbdir {}
	set dbfiles {}
	set realign 1
	set paired 1
	set adapterfile {}
	set conv_nextseq 0
	set skipincomplete 1
	set targetsfile {}
	set addtargets 0
	cg_options process_multicompar args {
		-dbdir {
			set dbdir $value
		}
		-split {
			set split $value
		}
		-dbfile {
			lappend dbfiles $value
		}
		-dbfiles {
			lappend dbfiles {*}$value
		}
		-skipincomplete {
			set skipincomplete $value
		}
		-targetsfile {
			set addtargets 1
			set targetsfile $value
		}
		-todo {
			set todo $value
		}
		-varfiles {
			set varfiles $value
		}
		-experiment {
			set experiment $value
		}
		-m - --maxopenfiles {
			set ::maxopenfiles [expr {$value - 4}]
		}
	} {destdir dbdir todo} 1 3
	set destdir [file_absolute $destdir]
	file mkdir $destdir/compar
	if {![info exists experiment]} {
		set experiment [file tail $destdir]
	}
	# check projectinfo
	projectinfo $destdir dbdir {split 0}
	set dbdir [dbdir $dbdir]
	set refseq [glob $dbdir/genome_*.ifas]
	# analysis info
	# -------------
	info_analysis_file $destdir/compar/info_analysis.tsv {} \
		{dbdir split dbfiles targetsfile ::maxopenfiles} \
		{genomecomb dbdir gnusort8 lz4 os} \
		command [list cg process_multicompar {*}$keepargs]

	# get samples todo
	# ----------------
	set samples {}
	if {[file exists $destdir/samples]} {
		set sampledir $destdir/samples
	} else {
		set sampledir $destdir
	}
	if {![info exists todo]} {
		set todo {}
		if {![info exists varfiles]} {
			set varfiles [jobglob ${sampledir}/*/var-*.tsv]
		} elseif {[string first * $varfiles] != -1} {
			set varfiles [jobglob $varfiles]
		}
		foreach file [lsort -dict $varfiles] {
			lappend todo [string range [file root [file tail [gzroot $file]]] 4 end]
		}
	}
	set todo [list_remdup $todo]
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
	putslog "Finding samples"
	if {[catch {cg select -n [gzfile compar/compar-$experiment.tsv]} done]} {set done {}}
	set compar_file compar/compar-$experiment.tsv.lz4
	set done [split $done \n]
	set stilltodo {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		if {![inlist $done $sample]} {
			putslog "Still todo: $sample"
			lappend stilltodo $sampledir/$name/var-$sample.tsv
		}
	}
	putslog "Samples: [llength $todo] todo, [llength $done] already done, [llength $stilltodo] to add"
	if {[llength $stilltodo]} {putslog "Samples to add: $stilltodo"}
	if {[llength $stilltodo] || $addtargets} {
		putslog "Starting multicompar"
		pmulticompar_job $compar_file $stilltodo 0 $split $targetsfile 0 $skipincomplete
	}
	# annotate multicompar
	# --------------------
	putslog "Starting annotation"
	cg_annotate_job $compar_file compar/annot_compar-$experiment.tsv.lz4 $dbdir {*}$dbfiles
	job indexannotcompar-$experiment \
	-deps compar/annot_compar-$experiment.tsv \
	-targets compar/annot_compar-$experiment.tsv.index/info.tsv -vars dbdir -code {
		cg index -colinfo $dep
	}
	#
	# multi sreg
	# ----------
	putslog "Starting multisreg"
	set regfiles {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		set file $sampledir/$name/sreg-$sample.tsv
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
	multireg_job compar/sreg-$experiment.tsv.lz4 $regfiles
	#
	# cgsv
	# ----
	putslog "Starting cgsv"
	set files [jobglob $sampledir/*/cgsv-*.tsv]
	if {[llength $files]} {
		set target compar/cgsv-${experiment}.tsv.lz4
		set names [list_regsub {.*/cgsv-(.*)\.tsv.*} $files {\1}]
		testmultitarget $target $names "$sampledir/\$name/cgsv-\$name.tsv"
		job cgsv_multicompar -optional 1 -deps $files -targets {$target} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -n $target.temp]
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
			cg lz4 $target.temp
			file rename -force $target.temp.lz4 $target
		}
		job cgsv_annotate -optional 1 \
		-deps {compar/cgsv-$experiment.tsv} \
		-targets {compar/annot_cgsv-$experiment.tsv.lz4} -vars {refseqdir build} -code {
			cg annotate $dep $target $refseqdir/$build
		}
		job cgsv_annotate_index -optional 1 \
		-deps {compar/annot_cgsv-$experiment.tsv.lz4} \
		-targets {compar/annot_cgsv-$experiment.tsv.index/info.tsv} -code {
			cg index $dep
		}
	}
	# cgcnv
	# ----
	putslog "Starting cgcnv"
	set files [jobglob $sampledir/*/cgcnv-*.tsv]
	if {[llength $files]} {
		set target compar/cgcnv-${experiment}.tsv.lz4
		set names [list_regsub {.*/cgcnv-(.*)\.tsv.*} $files {\1}]
		testmultitarget $target $names "$sampledir/\$name/cgcnv-\$name.tsv"
		job cgcnv_multicompar -optional 1 -deps $files -targets {compar/cgcnv-${experiment}.tsv.lz4} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -n $target.temp]
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
			cg lz4 $target.temp
			file rename -force $target.temp.lz4 $target
		}
		job cgcnv_annotate -optional 1 \
		-deps {compar/cgcnv-$experiment.tsv.lz4} \
		-targets {compar/annot_cgcnv-$experiment.tsv.lz4} -vars {refseqdir build} -code {
			cg annotate $dep $target $refseqdir/$build
		}
		job cgcnv_annotate_index -optional 1 \
		-deps {compar/annot_cgcnv-$experiment.tsv.lz4} \
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

