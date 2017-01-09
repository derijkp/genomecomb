proc process_multicompar_job {destdir experiment dbdir todo args} {
# putsvars destdir experiment dbdir todo args
	upvar job_logdir job_logdir
	set keeppwd [pwd]
	cd $destdir
	set skipincomplete 1
	set split 0
	set dbfiles {}
	set addtargets 0
	set targetsfile {}
	foreach {key value} $args {
		if {$key eq "-skipincomplete"} {
			set skipincomplete $value
		} elseif {$key eq "-split"} {
			set split $value
		} elseif {$key eq "-dbfiles"} {
			set dbfiles $value
		} elseif {$key eq "-targetsfile"} {
			set addtargets 1
			set targetsfile $value
		} else {
			lappend opts $key $value
		}
	}
	file mkdir compar
	#
	# multicompar
	putslog "Making multicompar"
	putslog "Finding samples"
	set compar_file compar/compar-$experiment.tsv
	if {[file exists samples]} {set sampledir samples/} else {set sampledir {}}
	if {[catch {cg select -n $compar_file} done]} {set done {}}
	set done [split $done \n]
	set stilltodo {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		if {![inlist $done $sample]} {
			putslog "Still todo: $sample"
			lappend stilltodo $sampledir$name/var-$sample.tsv
		}
	}
	if {[llength $stilltodo] || $addtargets} {
		putslog "Starting multicompar"
		pmulticompar_job $compar_file $stilltodo 0 $split $targetsfile 0 $skipincomplete
	}
	putslog "Starting annotation"
	job annotcompar-$experiment -deps [list $compar_file {*}$dbfiles] \
	-targets compar/annot_compar-$experiment.tsv -vars {dbdir dbfiles} -code {
		cg annotate $dep $target.temp $dbdir {*}$dbfiles
		file delete -force $target.temp.index
		file rename -force $target.temp $target
	}
	job indexannotcompar-$experiment \
	-deps compar/annot_compar-$experiment.tsv \
	-targets compar/annot_compar-$experiment.tsv.index/info.tsv -vars dbdir -code {
		cg index -colinfo $dep
	}
	#
	# multi sreg
	putslog "Starting multisreg"
	set regfiles {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		set file $sampledir$name/sreg-$sample.tsv
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
	multireg_job compar/sreg-$experiment.tsv $regfiles
	#
	# cgsv
	# ----
	putslog "Starting cgsv"
	set files [jobglob $sampledir*/cgsv-*.tsv]
	if {[llength $files]} {
		set target compar/cgsv-${experiment}.tsv
		set names [list_regsub {.*/cgsv-(.*)\.tsv.*} $files {\1}]
		testmultitarget $target $names "$sampledir\$name/cgsv-\$name.tsv"
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
			file rename -force $target.temp $target
		}
		job cgsv_annotate -optional 1 \
		-deps {compar/cgsv-$experiment.tsv} \
		-targets {compar/annot_cgsv-$experiment.tsv} -vars {refseqdir build} -code {
			cg annotate $dep $target.temp $refseqdir/$build
			file rename -force $target.temp $target
		}
		job cgsv_annotate_index -optional 1 \
		-deps {compar/annot_cgsv-$experiment.tsv} \
		-targets {compar/annot_cgsv-$experiment.tsv.index/info.tsv} -code {
			cg index $dep
		}
	}
	# cgcnv
	# ----
	putslog "Starting cgcnv"
	set files [jobglob $sampledir*/cgsv-*.tsv]
	if {[llength $files]} {
		set target compar/cgcnv-${experiment}.tsv
		set names [list_regsub {.*/cgsv-(.*)\.tsv.*} $files {\1}]
		testmultitarget $target $names "$sampledir\$name/cgcnv-\$name.tsv"
		job cgcnv_multicompar -optional 1 -deps $files -targets {compar/cgcnv-${experiment}.tsv} -code {
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
				puts "Multicgcnv already present: $done"
			}
			if {[llength $todo]} {
				cg svmulticompar $target.temp {*}$todo
			}
			file rename -force $target.temp $target
		}
		job cgcnv_annotate -optional 1 \
		-deps {compar/cgcnv-$experiment.tsv} \
		-targets {compar/annot_cgcnv-$experiment.tsv} -vars {refseqdir build} -code {
			cg annotate $dep $target.temp $refseqdir/$build
			file rename -force $target.temp $target
		}
		job cgcnv_annotate_index -optional 1 \
		-deps {compar/annot_cgcnv-$experiment.tsv} \
		-targets {compar/annot_cgcnv-$experiment.tsv.index/info.tsv} -code {
			cg index -colinfo $dep
		}
	}
	cd $keeppwd
}

proc process_multicompar {args} {
	set dbdir {}
	set dbfiles {}
	set realign 1
	set paired 1
	set adapterfile {}
	set conv_nextseq 0
	set skipincomplete 1
	set targetsfile {}
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
		-skipincomplete {
			set skipincomplete $value
		}
		-targetsfile {
			set targetsfile $value
		}
		-varfiles {
			set varfiles $value
		}
		-m - --maxopenfiles {
			set ::maxopenfiles [expr {$value - 4}]
		}
	} {destdir dbdir} 1 2
	set destdir [file_absolute $destdir]
	# check projectinfo
	projectinfo $destdir dbdir split
	set refseq [glob $dbdir/genome_*.ifas]
	set samples {}
	set experiment [file tail $destdir]
	if {[file exists $destdir/samples]} {
		set sampledir $destdir/samples
	} else {
		set sampledir $destdir
	}
	set todo {}
	if {![info exists varfiles]} {
		set varfiles [jobglob ${sampledir}/*/var-*.tsv]
	} elseif {[string first * $varfiles] != -1} {
		set varfiles [jobglob $varfiles]
	}
	foreach file [lsort -dict $varfiles] {
		lappend todo [string range [file root [file tail [gzroot $file]]] 4 end]
	}
	job_logdir $destdir/log_jobs
	set keeppwd [pwd]
	set todo [list_remdup $todo]
	process_multicompar_job $destdir $experiment $dbdir $todo \
		-skipincomplete $skipincomplete -split $split -dbfiles $dbfiles -targetsfile $targetsfile
}

proc cg_process_multicompar {args} {
	putslog "Running [list cg process_multicompar {*}$args]"
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_multicompar
	}
	process_multicompar {*}$args
	job_wait
}

