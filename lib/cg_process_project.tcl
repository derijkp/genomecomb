proc process_project_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg process_project {*}$args]
	unset -nocomplain split
	set preset {}
	set dbdir {}
	set dbfiles {}
	set organelles {}
	set minfastqreads {}
	set clip {}
	set removeskew {}
	set aligners bwa
	set ali_keepcomments {}
	set varcallers {}
	set isocallers {}
	set iso_joint {}
	set iso_match {}
	set counters {}
	set svcallers {}
	set methcallers {}
	set realign {}
	set cleanup {}
	set paired {}
	set adapterfile {}
	set conv_nextseq 0
	set targetfile {}
	set targetvarsfile {}
	set reports {}
	set samBQ {}
	set dt {}
	set removeduplicates {}
	set amplicons {}
	set extra_reports_mastr {}
	set jobsample 0
	set distrreg {}
	set threads {}
	set keepsams {}
	set keepfields *
	set datatype {}
	set aliformat {}
	set samplesdir samples
	set maxfastqdistr {}
	set hap_bam {}
	set depth_histo_max {}
	set reftranscripts {}
	set singlecell {}
	set sc_whitelist {}
	set sc_umisize {}
	set sc_barcodesize {}
	set sc_adaptorseq {}
	set sc_filters {}
	set sc_celltypers {}
	set sc_expectedcells {}
	set cellmarkerfile {}
	set tissue {}
	set optionsfile options.tsv
	set process_msamples 0
	cg_options process_project args {
		-preset {
			if {$value ne ""} {
				if {![command_exists preset_$value]} {
					error "preset $value does not exist, must be one of: [presets]"
				}
				foreach {var val} [preset_$value] {
					set $var $val
				}
				set preset $value
			}
		}
		-ori {
			set oridir $value
		}
		-dbdir - -refdir {
			set dbdir $value
#			set refseq [refseq $dbdir]
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
			set aligners $value
		}
		-ali_keepcomments {
			set ali_keepcomments $value
		}
		-singlecell {
			if {$value ni {ontr10x {}}} {error "Unknown value $value for -singlecell, must be one of: ontr10x (or empty)"}
			set singlecell $value
		}
		-sc_whitelist {
			set sc_whitelist $value
		}
		-sc_umisize {
			set sc_umisize $value
		}
		-sc_barcodesize {
			set sc_barcodesize $value
		}
		-sc_adaptorseq {
			set sc_adaptorseq $value
		}
		-sc_filters {
			set sc_filters $value
		}
		-sc_celltypers {
			set sc_celltypers $value
		}
		-sc_expectedcells {
			set sc_expectedcells $value
		}
		-cellmarkerfile {
			if {$value ne "" && ![file exists $value]} {error "cellmarkerfile $value does not exists"}
			set cellmarkerfile [file_absolute $value]
		}
		-tissue {
			set tissue $value
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
		-isocallers {
			set isocallers $value
		}
		-organelles {
			set organelles $value
		}
		-iso_joint {
			set iso_joint $value
		}
		-reftranscripts {
			set reftranscripts $value
		}
		-iso_match {
			set iso_match $value
		}
		-flair {
			if {[true $value]} {list_addnew isocallers flair}
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
		-process_msamples {
			set process_msamples $value
		}
		-*-* {
			set ::specialopt($key) $value
		}
	} {destdir dbdir} 1 2
	set destdir [file_absolute $destdir]
	set adapterfile [adapterfile $adapterfile]
	set experimentname [file tail $destdir]
	if {[file pathtype $optionsfile] ne "absolute"} {
		set optionsfile [file join $destdir $optionsfile]
	}
	unset -nocomplain optionsa
	if {[file exists $optionsfile]} {
		set f [gzopen $optionsfile]
		set header [tsv_open $f]
		if {$header ne "sample option value"} {error "wrong fields in optionsfile $optionsfile: must be: sample option value"}
		while {[gets $f line] != -1} {
			set sline [split $line \t]
			if {[llength $sline] != 3} {error "format error in optionsfile $optionsfile at line:$line"}
			foreach {sample option value} $sline break
			set optionsa($sample,$option) $value
			puts "sample options: $sample\t$option\t$value"
		}
		gzclose $f
	}
	# check projectinfo
	set dbdir [dbdir $dbdir]
	projectinfo $destdir dbdir {split 1}
	set ref [file tail $dbdir]
	if {$cellmarkerfile ne ""} {
		if {$sc_celltypers eq ""} {set sc_celltypers {scsorter sctype}}
	} elseif {$tissue ne ""} {
		if {$sc_celltypers eq ""} {set sc_celltypers {sctype}}
	}
	# logfile
	# -------
	job_logfile $destdir/process_project_[file tail $destdir] $destdir $cmdline \
		{*}[versions dbdir fastq-stats samtools gnusort8 zst os]
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
		if {$process_msamples} {
			foreach dir [jobglob $destdir/msamples/*] {
				set sample [file tail $dir]
				if {[regexp {[- ]} $sample]} {
					error "incompatible sample name $sample: sample names cannot contain spaces or dashes (-)"
				}
				lappend msamples $sample
			}
		}
	} else {
		set sampledir $destdir
		foreach dir [jobglob $sampledir/*/fastq $sampledir/*/ubam $sampledir/*/ori] {
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
	set keys {
		clip 
		singlecell sc_whitelist sc_umisize sc_barcodesize sc_adaptorseq 
		sc_filters sc_celltypers sc_expectedcells cellmarkerfile tissue 
		datatype aliformat aligners ali_keepcomments realign 
		varcallers svcallers methcallers counters reftranscripts isocallers 
		organelles hap_bam dbdir split paired maxfastqdistr adapterfile 
		reports samBQ cleanup removeduplicates amplicons threads distrreg 
		keepsams removeskew dt targetfile minfastqreads depth_histo_max 
	}
	foreach sample $samples {
		putslog "Processing sample $sample"
		set dir $sampledir/$sample
		set sampleargs [list \
			-preset [get optionsa($sample,preset) $preset] \
		]
		foreach key $keys {
			if {[info exists optionsa($sample,$key)]} {
				lappend sampleargs -$key $optionsa($sample,$key)
			} else {
				set value [get $key]
				if {$value ne ""} {
					if {$value ne ""} {
						lappend sampleargs -$key $value
					}
				}
			}
		}
		if {!$jobsample} {
			process_sample_job -todoVar todo {*}$sampleargs $dir
		} else {
			# find deps and targets by running the process_sample_job with job_getinfo set to 1
			job_getinfo 1
			set verbose [logverbose]
			set ::deps {} ; set ::targets {}
			process_sample_job -todoVar todo {*}$sampleargs $dir
			foreach {deps targets} [job_getinfo 0] break
			if {$targetfile ne ""} {lappend deps $targetfile}
			logverbose $verbose
			# run the actual job with deps and targets found
			job process_sample-$sample -deps $deps -targets $targets -vars {
				sampleargs dir
			} -code {
				cg process_sample -stack 1 -v 2 {*}$sampleargs $dir >@ stdout 2>@ stderr
			}
		}
	}
	# run msamples (not supporting jobsample for these)
	foreach sample $msamples {
		putslog "Processing msample $sample"
		set dir $destdir/msamples/$sample
		set sampleargs [list \
			-preset [get optionsa($sample,preset) $preset] \
		]
		foreach key $keys {
			set value [get $key]
			if {$value ne ""} {
				if {[info exists optionsa($sample,$key)]} {
					lappend sampleargs -$key $optionsa($sample,$key)
				} elseif {$value ne ""} {
					lappend sampleargs -$key $value
				}
			}
		}
		process_sample_job {*}$sampleargs $dir
	}
	set_job_logdir $destdir/log_jobs
	set todo(var) [list_remdup $todo(var)]
	set todo(sv) [list_remdup $todo(sv)]
	set todo(meth) [list_remdup $todo(meth)]
	set todo(reports) [list_remdup $todo(reports)]
	if {![llength $reports]} {set todo(reports) {}}
	if {$cleanup eq ""} {set cleanup 1}
	process_multicompar_job \
		-experiment $experiment \
		-skipincomplete 1 -targetvarsfile $targetvarsfile \
		-isocallers $isocallers \
		-iso_match $iso_match \
		-sc_celltypers $sc_celltypers \
		-threads $threads \
		-distrreg $distrreg \
		-keepfields $keepfields \
		-split $split -dbfiles $dbfiles \
		-cleanup $cleanup \
		$destdir $dbdir
	if {$extra_reports_mastr ne ""} {
		make_alternative_compar_job $experiment $destdir $extra_reports_mastr
		set histofiles [jobglob $destdir/samples/*/reports/*.histo]
		generate_coverage_report_job $experiment $amplicons $histofiles $destdir
		generate_html_report_job $experiment $destdir
		analysis_complete_job $experiment $destdir $extra_reports_mastr
	}
	if {[llength $iso_joint]} {
		iso_joint_job \
			-iso_joint $iso_joint \
			-iso_match $iso_match \
			-threads $threads \
			-distrreg $distrreg \
			-dbdir $dbdir \
			-cleanup $cleanup \
			$destdir
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
