proc process_project_job {args} {
	set dbdir {}
	set dbfiles {}
	set realign 1
	set paired 1
	set adapterfile {}
	set conv_nextseq 0
	set reports all
	cg_options process_project_job args {
		-ori {
			set oridir $value
		}
		-realign {
			set realign $value
		}
		-dbdir - -refdir {
			set dbdir $value
		}
		-split {
			set split $value
		}
		-dbfile {
			lappend dbfiles $value
		}
		-paired {
			set paired $value
		}
		-adapterfile {
			set adapterfile $value
		}
		-conv_nextseq {
			set conv_nextseq $value
		}
		-reports {
			set reports $value
		}
	} 1 2
	set len [llength $args]
	if {$len == 1} {
		set destdir [lindex $args 0]
	} elseif {$len == 2} {
		foreach {destdir dbdir} $args break
	} else {
		errorformat process_illumina
	}
	set destdir [file_absolute $destdir]
	set dbdir [file_absolute $dbdir]
	# check projectinfo
	projectinfo $destdir dbdir {split 1}
	# start
	##in case of nextseq500 data - generate fastqs & distribute data
	if {$conv_nextseq} {
		set rundir [glob $destdir/*NS500*]
		cg_bcl2fastq $rundir fastq 4 4 4 4
		cg_process_conv_illnextseq fastq $destdir
	}
	set samples {}
	set experiment [file tail $destdir]
	if {[file exists $destdir/samples]} {
		set sampledir $destdir/samples
	} else {
		set sampledir $destdir
	}
	unset -nocomplain a
	foreach dir [jobglob $sampledir/*/fastq $sampledir/*/ori] {
		set a([file tail [file dir $dir]]) 1
	}
	set samples [ssort -natural [array names a]]
	set keeppwd [pwd]
	cd $destdir
	job_logdir $destdir/log_jobs
	set todo {}
	set reportstodo {}
	foreach sample $samples {
		putslog "Processing sample $sample"
		set dir $sampledir/$sample
		process_sample_job -todoVar todo -reportstodoVar reportstodo \
			-realign $realign -dbdir $dbdir -split $split -paired $paired \
			-adapterfile $adapterfile -reports $reports \
			$dir

	}
	job_logdir $destdir/log_jobs
	cd $destdir
	set todo [list_remdup $todo]
	process_multicompar_job $destdir $experiment $dbdir $todo -skipincomplete 1 -split $split -dbfiles $dbfiles
	if {[llength $reports]} {
		proces_reportscombine_job $destdir $reportstodo
	}
	cd $keeppwd
	list $todo $
}

proc cg_process_project {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_project
	}
	process_project_job {*}$args
	job_wait
}
