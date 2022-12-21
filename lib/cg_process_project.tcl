proc process_project_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg process_project {*}$args]
	unset -nocomplain split
	set dbdir {}
	set dbfiles {}
	set minfastqreads 0
	set clip 1
	set removeskew {}
	set aligner bwa
	set varcallers {gatkh strelka}
	set counters {}
	set svcallers {}
	set methcallers {}
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
	set extra_reports_mastr {}
	set jobsample 0
	set distrreg 0
	set threads 2
	set keepsams 0
	set keepfields *
	set datatype {}
	set aliformat bam
	set samplesdir samples
	set maxfastqdistr {}
	set hap_bam 0
	set depth_histo_max 1000
	set flair 0
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
		-maxfastqdistr {
			set maxfastqdistr $value
		}
		-clip {
			set clip $value
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
		-a - -aligner - -aligners {
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
		-svcallers {
			set svcallers $value
		}
		-methcallers {
			set methcallers $value
		}
		-counters {
			set counters $value
		}
		-flair {
			set flair $value
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
		-keepfields {
			set keepfields $value
		}
		-r - -reports {
			set reports $value
		}
		-c - -cleanup {
			set cleanup $value
		}
		-m - -maxopenfiles {
			maxopenfiles [expr {$value - 4}]
		}
		-threads {
			set threads $value
		}
		-distrreg {
			set distrreg [distrreg_checkvalue $value]
		}
		-conv_nextseq {
			set conv_nextseq $value
		}
		-extra_reports_mastr {
			if {$value eq "1"} {set value gatk}
			set extra_reports_mastr $value
		}
		-jobsample {
			if {$value ni {0 1}} {error "-jobsample must be 0 or 1"}
			set jobsample $value
		}
		-keepsams {
			set keepsams $value
		}
		-datatype {
			set datatype $value
		}
		-aliformat {
			set aliformat $value
		}
		-samplesdir {
			set samplesdir $value
		}
		-hap_bam {
			set hap_bam $value
		}
		-depth_histo_max {
			set depth_histo_max $value
		}
		-extraopts {
			foreach {k v} $value {
				set ::cgextraopts($k) $v
			}
		}
		-*-* {
			set ::specialopt($key) $value
		}
	} {destdir dbdir} 1 2
	set destdir [file_absolute $destdir]
	set adapterfile [adapterfile $adapterfile]
	set experimentname [file tail $destdir]
	# check projectinfo
	set dbdir [dbdir $dbdir]
	projectinfo $destdir dbdir {split 1}
	set ref [file tail $dbdir]
	# logfile
	# -------
	job_logfile $destdir/process_project_[file tail $destdir] $destdir $cmdline \
		{*}[versions dbdir fastqc fastq-stats fastq-mcf bwa samtools gatk gatk3 picard java gnusort8 zst os]
	# amplicons settings
	# ------------------
	if {$amplicons ne ""} {
		if {$removeduplicates eq ""} {set removeduplicates 0}
		if {$removeskew eq ""} {set removeskew 0}
		if {$dt eq ""} {set dt NONE}
		list_addnew dbfiles $amplicons
		if {$targetfile eq ""} {
			set targetfile $destdir/reg_${ref}_targets.tsv.zst
			job reports_amplicons2targetfile -deps {$amplicons} -targets {$targetfile} -vars {sample dbdir ref} -code {
				cg regcollapse $dep | cg zst > $target
			}
		}
	} else {
		if {$removeduplicates eq ""} {set removeduplicates 1}
		if {$removeskew eq ""} {set removeskew 1}
	}
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
	set msamples {}
	if {[file exists $destdir/$samplesdir]} {
		set sampledir $destdir/$samplesdir
		foreach dir [jobglob $sampledir/*] {
			set sample [file tail $dir]
			if {[regexp {[- ]} $sample]} {
				error "incompatible sample name $sample: sample names cannot contain spaces or dashes (-)"
			}
			set a($sample) 1
		}
		foreach dir [jobglob $destdir/msamples/*] {
			set sample [file tail $dir]
			if {[regexp {[- ]} $sample]} {
				error "incompatible sample name $sample: sample names cannot contain spaces or dashes (-)"
			}
			lappend msamples $sample
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
	set msamples [bsort $msamples]
	set samples [bsort [array names a]]
	set poss [list_find -regexp $samples -]
	if {[llength $poss]} {
		error "- is not allowed in sample names. The following sample name(s) have a -: [list_sub $samples $poss]"
	}
	set_job_logdir $destdir/log_jobs
	set todo(var) {}
	set todo(sv) {}
	set todo(meth) {}
	set todo(reports) {}
	foreach sample $samples {
		putslog "Processing sample $sample"
		set dir $sampledir/$sample
		if {!$jobsample} {
			process_sample_job -todoVar todo -clip $clip -datatype $datatype -aliformat $aliformat \
				-aligner $aligner -realign $realign \
				-varcallers $varcallers -svcallers $svcallers -methcallers $methcallers \
				-counters $counters \
				-hap_bam $hap_bam \
				-dbdir $dbdir -split $split -paired $paired --maxfastqdistr $maxfastqdistr \
				-adapterfile $adapterfile -reports $reports -samBQ $samBQ -cleanup $cleanup \
				-removeduplicates $removeduplicates -amplicons $amplicons \
				-threads $threads -distrreg $distrreg -keepsams $keepsams \
				-removeskew $removeskew -dt $dt -targetfile $targetfile -minfastqreads $minfastqreads \
				$dir
		} else {
			# find deps and targets by running the process_sample_job with job_getinfo set to 1
			job_getinfo 1
			set verbose [logverbose]
			set ::deps {} ; set ::targets {}
			process_sample_job -todoVar todo -clip $clip -datatype $datatype \
				-aligner $aligner -realign $realign \
				-varcallers $varcallers -svcallers $svcallers -methcallers $methcallers \
				-counters $counters \
				-hap_bam $hap_bam \
				-dbdir $dbdir -split $split -paired $paired -keepsams $keepsams --maxfastqdistr $maxfastqdistr \
				-adapterfile $adapterfile -reports $reports -samBQ $samBQ -cleanup $cleanup \
				-removeduplicates $removeduplicates -amplicons $amplicons \
				-removeskew $removeskew -dt $dt -targetfile $targetfile -minfastqreads $minfastqreads\
				$dir
			foreach {deps targets} [job_getinfo 0] break
			logverbose $verbose
			# run the actual job with deps and targets found
			job process_sample-$sample -deps $deps -targets $targets -vars {
				clip aligner realign varcallers svcallers methcallers dbdir split paired
				adapterfile reports samBQ cleanup removeduplicates amplicons
				removeskew dt targetfile minfastqreads dir keepsams datatype maxfastqdistr
				counters
			} -code {
				cg process_sample -stack 1 -v 2 -clip $clip -datatype $datatype \
					-aligner $aligner -realign $realign \
					-varcallers $varcallers -svcallers $svcallers -methcallers $methcallers \
					-counters $counters \
					-dbdir $dbdir -split $split -paired $paired -keepsams $keepsams --maxfastqdistr $maxfastqdistr \
					-adapterfile $adapterfile -reports $reports -samBQ $samBQ -cleanup $cleanup \
					-removeduplicates $removeduplicates -amplicons $amplicons \
					-removeskew $removeskew -dt $dt -targetfile $targetfile -minfastqreads $minfastqreads\
					$dir >@ stdout 2>@ stderr
			}
		}
	}
	# run msamples (not supporting jobsample for these)
	foreach sample $msamples {
		putslog "Processing msample $sample"
		set dir $destdir/msamples/$sample
		process_sample_job -clip $clip -datatype $datatype -aliformat $aliformat \
			-aligner $aligner -realign $realign \
			-varcallers $varcallers -svcallers $svcallers -methcallers $methcallers \
			-counters $counters \
			-hap_bam $hap_bam \
			-dbdir $dbdir -split $split -paired $paired --maxfastqdistr $maxfastqdistr \
			-adapterfile $adapterfile -reports $reports -samBQ $samBQ -cleanup $cleanup \
			-removeduplicates $removeduplicates -amplicons $amplicons \
			-threads $threads -distrreg $distrreg -keepsams $keepsams \
			-removeskew $removeskew -dt $dt -targetfile $targetfile -minfastqreads $minfastqreads \
			-depth_histo_max $depth_histo_max \
			$dir
	}
	set_job_logdir $destdir/log_jobs
	set todo(var) [list_remdup $todo(var)]
	set todo(sv) [list_remdup $todo(sv)]
	set todo(meth) [list_remdup $todo(meth)]
	set todo(reports) [list_remdup $todo(reports)]
	if {![llength $reports]} {set todo(reports) {}}
	process_multicompar_job -experiment $experiment \
		-skipincomplete 1 -targetvarsfile $targetvarsfile \
		-varfiles $todo(var) -svfiles $todo(sv) -methfiles $todo(meth) \
		-counters $counters \
		-threads $threads -distrreg $distrreg \
		-keepfields $keepfields \
		-split $split -dbfiles $dbfiles -cleanup $cleanup \
		-reports $todo(reports) \
		$destdir $dbdir
	if {$extra_reports_mastr ne ""} {
		make_alternative_compar_job $experiment $destdir $extra_reports_mastr
		set histofiles {}
		foreach dir $todo(reports) {
			set list [jobglob $dir/*.histo]
			if {[llength $list]} {lappend histofiles {*}$list}
		}
		generate_coverage_report_job $experiment $amplicons $histofiles $destdir
		generate_html_report_job $experiment $destdir
		analysis_complete_job $experiment $destdir $extra_reports_mastr
	}
	if {$flair} {
		flair_job -refseq [refseq $dbdir] $destdir
	}
	list $todo(var)
}

proc cg_process_project {args} {
	set args [job_init {*}$args]
	putslog pwd:\ [pwd]
	putslog cmd:\ [list cg process_project {*}$args]
	putslog jobargs:\ [job_curargs]
	if {[llength $args] < 1} {
		errorformat process_project
	}
	process_project_job {*}$args
	job_wait
}
