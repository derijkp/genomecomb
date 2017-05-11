proc process_project_job {args} {
	set dbdir {}
	set dbfiles {}
	set aligner bwa
	set varcallers {gatk sam}
	set realign 1
	set cleanup 1
	set paired 1
	set adapterfile {}
	set conv_nextseq 0
	set reports all
	set samBQ 0
	cg_options process_project_job args {
		-ori {
			set oridir $value
		}
		-dbdir - -refdir {
			set dbdir $value
		}
		-a - -aligner {
			set aligner $value
		}
		-realign {
			set realign $value
		}
		-v - -varcallers {
			set varcallers $value
		}
		-s - -split {
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
				if {![file exists $file]} {error "dbfile $v does not exists"}
				lappend dbfiles [file_absolute $file]
			}
		}
		-p - -paired {
			set paired $value
		}
		-adapterfile {
			if {$value ne "" && ![file exists $value]} {error "adapterfile $value does not exists"}
			set adapterfile $value
		}
		-conv_nextseq {
			set conv_nextseq $value
		}
		-r - -reports {
			set reports $value
		}
		-c - -cleanup {
			set cleanup $value
		}
		-m - -maxopenfiles {
			set ::maxopenfiles [expr {$value - 4}]
		}
		-samBQ {
			set samBQ $value
		}
	} {destdir dbdir} 1 2
	set destdir [file_absolute $destdir]
	set dbdir [file_absolute $dbdir]
	# check projectinfo
	projectinfo $destdir dbdir {split 1}
	set dbdir [dbdir $dbdir]
	# logfile
	# -------
	set cmdline [list cg process_project]
	foreach option {
		ori dbdir refdir a aligner realign v varcallers split dbfile dbfiles paired adapterfile conv_nextseq reports cleanup maxopenfiles samBQ
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline $destdir $dbdir
	job_logfile $destdir/process_project_[file tail $destdir] $destdir $cmdline \
		{*}[versions dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk picard java gnusort8 lz4 os]
	# start
	# -----
	##in case of nextseq500 data - generate fastqs & distribute data
	if {$conv_nextseq} {
		set rundir [glob $destdir/*NS500*]
		cg_bcl2fastq $rundir fastq 4 4 4 4
		cg_process_conv_illnextseq fastq $destdir
	}
	set samples {}
	set experiment [file tail $destdir]
	unset -nocomplain a
	if {[file exists $destdir/samples]} {
		set sampledir $destdir/samples
		foreach dir [jobglob $sampledir/*] {
			set a([file tail $dir]) 1
		}
	} else {
		set sampledir $destdir
		foreach dir [jobglob $sampledir/*/fastq $sampledir/*/ori] {
			set a([file tail [file dir $dir]]) 1
		}
	}
	set samples [ssort -natural [array names a]]
	set poss [list_find -regexp $samples -]
	if {[llength $poss]} {
		error "- is not allowed in sample names. The following sample name(s) have a -: [list_sub $samples $poss]"
	}
	job_logdir $destdir/log_jobs
	set todo {}
	set reportstodo {}
	foreach sample $samples {
		putslog "Processing sample $sample"
		set dir $sampledir/$sample
		process_sample_job -todoVar todo -reportstodoVar reportstodo \
			-aligner $aligner -realign $realign --varcallers $varcallers \
			-dbdir $dbdir -split $split -paired $paired \
			-adapterfile $adapterfile -reports $reports -samBQ $samBQ -cleanup $cleanup \
			$dir
	}
	job_logdir $destdir/log_jobs
	set todo [list_remdup $todo]
	process_multicompar_job -experiment $experiment -skipincomplete 1 \
		-split $split -dbfiles $dbfiles -cleanup $cleanup $destdir $dbdir $todo
	if {[llength $reports]} {
		proces_reportscombine_job $destdir $reportstodo
	}
	list $todo
}

proc cg_process_project {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_project
	}
	process_project_job {*}$args
	job_wait
}
