proc process_project_job {args} {
	set dbdir {}
	set dbfiles {}
	set minfastqreads 0
	set removeskew {}
	set aligner bwa
	set varcallers {gatk sam}
	set realign 1
	set cleanup 1
	set paired 1
	set adapterfile {}
	set conv_nextseq 0
	set targetfile {}
	set targetvarsfile {}
	set reports basic
	set samBQ 0
	set dt {}
	set removeduplicates {}
	set amplicons {}
	set extra_reports_mastr 0
	set jobsample 0
	cg_options process_project args {
		-ori {
			set oridir $value
		}
		-dbdir - -refdir {
			set dbdir $value
		}
		-minfastqreads {
			set minfastqreads $value
		}
		-p - -paired {
			set paired $value
		}
		-adapterfile {
			if {$value ne "" && ![file exists $value]} {error "adapterfile $value does not exists"}
			set adapterfile $value
		}
		-removeskew {
			set removeskew $value
		}
		-a - -aligner {
			set aligner $value
		}
		-realign {
			set realign $value
		}
		-removeduplicates {
			set removeduplicates $value
		}
		-amplicons {
			set amplicons [file_absolute $value]
			if {$value ne "" && ![jobfileexists $amplicons]} {error "amplicons file $amplicons does not exists"}
		}
		-v - -varcallers {
			set varcallers $value
		}
		-s - -split {
			set split $value
		}
		-dt - -downsampling_type {
			if {$value ni "{} NONE ALL_READS BY_SAMPLE"} {error "-dt must be one of: NONE ALL_READS BY_SAMPLE"}
			set dt $value
		}
		-samBQ {
			set samBQ $value
		}
		-dbfile {
			set file [gzfile $value]
			if {$value ne "" && ![file exists $file]} {error "dbfile $value does not exists"}
			lappend dbfiles [file_absolute $file]
		}
		-dbfiles {
			foreach v $value {
				set file [gzfile $v]
				if {![file exists $file]} {error "dbfile $v does not exists"}
				lappend dbfiles [file_absolute $file]
			}
		}
		-targetfile {
			set targetfile [file_absolute $value]
			if {$value ne "" && ![jobfileexists $targetfile]} {error "target file $targetfile does not exists"}
		}
		-targetvarsfile {
			set targetvarsfile [file_absolute $value]
			if {$targetvarsfile ne "" && ![jobfileexists $value]} {error "targetvarsfile $targetvarsfile does not exists"}
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
		-conv_nextseq {
			set conv_nextseq $value
		}
		-extra_reports_mastr {
			set extra_reports_mastr $value
		}
		-jobsample {
			if {$value ni {0 1}} {error "-jobsample must be 0 or 1"}
			set jobsample $value
		}
	} {destdir dbdir} 1 2
	set destdir [file_absolute $destdir]
	set dbdir [file_absolute $dbdir]
	set adapterfile [adapterfile $adapterfile]
	set experimentname [file tail $destdir]
	if {$amplicons ne ""} {
		if {$removeduplicates eq ""} {set removeduplicates 0}
		if {$removeskew eq ""} {set removeskew 0}
		if {$dt eq ""} {set dt NONE}
		list_addnew dbfiles $amplicons
	} else {
		if {$removeduplicates eq ""} {set removeduplicates 1}
		if {$removeskew eq ""} {set removeskew 1}
		if {$dt eq ""} {set dt BY_SAMPLE}
	}
	# check projectinfo
	projectinfo $destdir dbdir {split 1}
	set dbdir [dbdir $dbdir]
	# logfile
	# -------
	set cmdline [list cg process_project]
	foreach option {
		ori dbdir refdir aligner realign varcallers split dbfile dbfiles paired adapterfile conv_nextseq reports cleanup jobsample maxopenfiles samBQ amplicons extra_reports_mastr
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline $destdir $dbdir
	job_logfile $destdir/process_project_[file tail $destdir] $destdir $cmdline \
		{*}[versions dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk gatk3 biobambam picard java gnusort8 lz4 os]
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
			set sample [file tail $dir]
			if {[regexp {[- ]} $sample]} {
				error "incompatible sample name $sample: sample names cannot contain spaces or dashes (-)"
			}
			set a($sample) 1
		}
	} else {
		set sampledir $destdir
		foreach dir [jobglob $sampledir/*/fastq $sampledir/*/ori] {
			set sample [file tail [file dir $dir]]
			if {[regexp {[- ]} $sample]} {
				error "incompatible sample name $sample: sample names cannot contain spaces or dashes (-)"
			}
			set a($sample) 1
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
		if {!$jobsample} {
			process_sample_job -todoVar todo -reportstodoVar reportstodo \
				-aligner $aligner -realign $realign --varcallers $varcallers \
				-dbdir $dbdir -split $split -paired $paired \
				-adapterfile $adapterfile -reports $reports -samBQ $samBQ -cleanup $cleanup \
				-removeduplicates $removeduplicates -amplicons $amplicons \
				-removeskew $removeskew -dt $dt -targetfile $targetfile -minfastqreads $minfastqreads \
				$dir
		} else {
			# find deps and targets by running the process_sample_job with job_getinfo set to 1
			job_getinfo 1
			set verbose [logverbose]
			set ::deps {} ; set ::targets {}
			process_sample_job -todoVar todo -reportstodoVar reportstodo \
				-aligner $aligner -realign $realign --varcallers $varcallers \
				-dbdir $dbdir -split $split -paired $paired \
				-adapterfile $adapterfile -reports $reports -samBQ $samBQ -cleanup $cleanup \
				-removeduplicates $removeduplicates -amplicons $amplicons \
				-removeskew $removeskew -dt $dt -targetfile $targetfile -minfastqreads $minfastqreads\
				$dir
			foreach {deps targets} [job_getinfo 0] break
			logverbose $verbose
			# run the actual job with deps and targets found
			job process_sample-$sample \
				-deps $deps -targets $targets -vars {
				aligner realign varcallers dbdir split paired
				adapterfile reports samBQ cleanup  removeduplicates amplicons
				removeskew dt targetfile minfastqreads dir
			} -code {
				cg process_sample -stack 1 -v 2 \
					-aligner $aligner -realign $realign --varcallers $varcallers \
					-dbdir $dbdir -split $split -paired $paired \
					-adapterfile $adapterfile -reports $reports -samBQ $samBQ -cleanup $cleanup \
					-removeduplicates $removeduplicates -amplicons $amplicons \
					-removeskew $removeskew -dt $dt -targetfile $targetfile -minfastqreads $minfastqreads\
					$dir >@ stdout 2>@ stderr
			}
		}
	}
	job_logdir $destdir/log_jobs
	set todo [list_remdup $todo]
	set reportstodo [list_remdup $reportstodo]
	process_multicompar_job -experiment $experiment \
		-skipincomplete 1 -targetvarsfile $targetvarsfile \
		-split $split -dbfiles $dbfiles -cleanup $cleanup $destdir $dbdir $todo
	if {[llength $reports]} {
		proces_reportscombine_job $destdir/reports {*}$reportstodo
		if {[jobfileexists $destdir/reports/report_hsmetrics-${experimentname}.tsv]} {
			mklink $destdir/reports/report_hsmetrics-${experimentname}.tsv $destdir/${experimentname}_hsmetrics_report.tsv
		}
	}
	if {$extra_reports_mastr} {
		make_alternative_compar_job $experiment $destdir
		set histofiles {}
		foreach dir $reportstodo {
			set list [jobglob $dir/*.histo]
			if {[llength $list]} {lappend histofiles {*}$list}
		}
		generate_coverage_report_job $experiment $amplicons $histofiles $destdir
		generate_html_report_job $experiment $destdir
		analysis_complete_job $experiment $destdir
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
