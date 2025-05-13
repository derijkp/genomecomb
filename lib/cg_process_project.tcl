proc code_empty {value} {
	if {$value eq ""} {
		return -
	} else {
		return $value
	}
}

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
	set aligners {}
	set ali_keepcomments {}
	set varcallers {}
	set isocallers {}
	set iso_joint {}
	set iso_joint_min 2
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
	set addumis 0
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
	set samplesheet {}
	set validate 0
	set extraannot {}
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
		-dbdir - -refdir - -refseq {
			# set dbdir $value
			set refseq [refseq $value]
			set dbdir [refdir $refseq]
		}
		-minfastqreads {
			set minfastqreads $value
		}
		-maxfastqdistr {
			set maxfastqdistr [code_empty $value]
		}
		-clip {
			set clip $value
		}
		-p - -paired {
			set paired $value
		}
		-adapterfile {
			if {$value ne "" && ![file exists $value]} {error "adapterfile $value does not exists"}
			set adapterfile [code_empty $value]
		}
		-removeskew {
			set removeskew [code_empty $value]
		}
		-a - -aligner - -aligners {
			set aligners $value
		}
		-ali_keepcomments {
			set ali_keepcomments [code_empty $value]
		}
		-singlecell {
			if {$value ni {ontr10x {}}} {error "Unknown value $value for -singlecell, must be one of: ontr10x (or empty)"}
			set singlecell [code_empty $value]
		}
		-addumis {
			if {$value ni {0 1}} {error "Unknown value $value for -addumis, must be either 1 or 0"}
			set addumis [code_empty $value]
		}
		-sc_whitelist {
			set sc_whitelist [code_empty $value]
		}
		-sc_umisize {
			set sc_umisize [code_empty $value]
		}
		-sc_barcodesize {
			set sc_barcodesize [code_empty $value]
		}
		-sc_adaptorseq {
			set sc_adaptorseq [code_empty $value]
		}
		-sc_filters {
			set sc_filters [code_empty $value]
		}
		-sc_celltypers {
			set sc_celltypers [code_empty $value]
		}
		-sc_expectedcells {
			set sc_expectedcells [code_empty $value]
		}
		-cellmarkerfile {
			if {$value ne "" && ![file exists $value]} {error "cellmarkerfile $value does not exists"}
			set cellmarkerfile [code_empty [file_absolute $value]]
		}
		-tissue {
			set tissue [code_empty $value]
		}
		-realign {
			set realign $value
		}
		-removeduplicates {
			set removeduplicates [code_empty $value]
		}
		-amplicons {
			if {$value ne "-" && ![jobfileexists $value]} {error "amplicons file $amplicons does not exists"}
			set amplicons [code_empty [file_absolute $value]]
		}
		-v - -varcallers {
			set varcallers [code_empty $value]
		}
		-svcallers {
			set svcallers [code_empty $value]
		}
		-methcallers {
			set methcallers [code_empty $value]
		}
		-counters {
			set counters [code_empty $value]
		}
		-isocallers {
			set isocallers [code_empty $value]
		}
		-organelles {
			set organelles [code_empty $value]
		}
		-iso_joint {
			set iso_joint $value
		}
		-iso_joint_min {
			set iso_joint_min $value
		}
		-extraannot {
			if {[llength [list_remove $value AnnotSV]]} {
				error "Unknown value(s) given for -extraannot: [list_remove $value AnnotSV]"
			}
			set extraannot $value
		}
		-reftranscripts {
			set reftranscripts [code_empty $value]
		}
		-iso_match {
			set iso_match $value
		}
		-flair {
			if {[true $value]} {
				list_addnew isocallers flair
				set isocallers [list_remove $isocallers -]
			}
		}
		-s - -split {
			set split $value
		}
		-dt - -downsampling_type {
			if {$value ni "{} NONE ALL_READS BY_SAMPLE"} {error "-dt must be one of: NONE ALL_READS BY_SAMPLE"}
			set dt [code_empty $value]
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
			if {$value ne "" && ![jobfileexists $value]} {error "target file $value does not exists"}
			set targetfile [code_empty [file_absolute $value]]
		}
		-targetvarsfile {
			if {$value ne "" && ![jobfileexists $value]} {error "targetvarsfile $value does not exists"}
			set targetvarsfile [file_absolute $value]
		}
		-keepfields {
			set keepfields $value
		}
		-r - -reports {
			set reports [code_empty $value]
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
			if {![isint $value]} {error "-jobsample must be an integer"}
			set jobsample $value
		}
		-keepsams {
			set keepsams $value
		}
		-datatype {
			set datatype [code_empty $value]
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
		-samplesheet {
			set samplesheet $value
		}
		-validate {
			set validate $value
			if {[true $validate]} {
				job_getinfo 1
				# logverbose 0
			}
		}
		-*-* {
			set ::specialopt($key) $value
		}
	} {destdir dbdir} 1 2
	set destdir [file_absolute $destdir]
	mkdir $destdir
	if {!$validate && $samplesheet ne ""} {
		cg make_project $destdir $samplesheet
	}
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
	if {!$validate} {
		job_logfile $destdir/process_project_[file tail $destdir] $destdir $cmdline \
			{*}[versions dbdir fastq-stats samtools gnusort8 zst os]
	}
	# amplicons settings
	# ------------------
	if {$amplicons ne ""} {
		list_addnew dbfiles $amplicons
		if {$targetfile eq ""} {
			set targetfile $destdir/reg_${ref}_targets.tsv.zst
			job reports_amplicons2targetfile -deps {$amplicons} -targets {$targetfile} -vars {sample dbdir ref} -code {
				cg regcollapse $dep | cg zst > $target
			}
		}
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
		singlecell addumis sc_whitelist sc_umisize sc_barcodesize sc_adaptorseq 
		sc_filters sc_celltypers sc_expectedcells cellmarkerfile tissue 
		datatype aliformat aligners ali_keepcomments realign 
		varcallers svcallers methcallers counters reftranscripts isocallers 
		organelles hap_bam dbdir split paired maxfastqdistr adapterfile 
		reports samBQ cleanup removeduplicates amplicons threads distrreg 
		keepsams removeskew dt targetfile minfastqreads depth_histo_max validate
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
		if {!$jobsample || $validate} {
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
			job process_sample-$sample -cores $jobsample -deps $deps -targets $targets -vars {
				sampleargs dir jobsample
			} -code {
				if {$jobsample > 1} {
					cg process_sample -stack 1 -v 2 -shadowdir [get ::env(SHADOWDIR) {}] -d $jobsample {*}$sampleargs $dir >@ stdout 2>@ stderr
				} else {
					cg process_sample -stack 1 -v 2 {*}$sampleargs $dir >@ stdout 2>@ stderr
				}
			}
		}
	}
	if {$validate} return
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
		-extraannot $extraannot \
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
		set reftranscripts [ref_tsvtranscripts $dbdir]
		iso_joint_job \
			-reftranscripts $reftranscripts \
			-iso_joint $iso_joint \
			-iso_joint_min $iso_joint_min \
			-iso_match $iso_match \
			-addumis $addumis \
			-threads $threads \
			-distrreg $distrreg \
			-dbdir $dbdir \
			-cleanup $cleanup \
			$destdir
	}
	list $todo(var)
}

proc cg_process_project {args} {
	# validation run
	set args [args_remove $args {-validate}]
	set keeppwd [pwd]
	set valargs [job_init {*}$args]
	putslog "validating process_sample $valargs"
	process_project_job -validate 1 {*}$valargs
	putslog "validation done. Start for real"
	# real run
	cd $keeppwd
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
