proc process_project_job {args} {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg process_project {*}$args]
	unset -nocomplain split
	set preset {}
	set dbdir {}
	set dbfiles {}
	set organelles {}
	set minfastqreads 0
	set clip 1
	set removeskew {}
	set aligners bwa
	set ali_keepcomments {}
	set varcallers {gatkh strelka}
	set isocallers {}
	set iso_joint {}
	set iso_match {}
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
	set reftranscripts {}
	set singlecell {}
	set singlecell_whitelist {}
	set singlecell_umisize 12
	set sc_filters {}
	set sc_celltypers {}
	set sc_expectedcells {}
	set cellmarkerfile {}
	set tissue {}
	set optionsfile options.tsv
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
			set ali_keepcomments [true $value]
		}
		-singlecell {
			if {$value ni {ontr10x {}}} {error "Unknown value $value for -singlecell, must be one of: ontr10x (or empty)"}
			set singlecell $value
		}
		-sc_whitelist {
			set singlecell_whitelist $value
		}
		-sc_umisize {
			set singlecell_umisize $value
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
	foreach sample $samples {
		putslog "Processing sample $sample"
		set dir $sampledir/$sample
		set sampleargs [list \
			-preset [get optionsa($sample,preset) $preset] \
			-clip [get optionsa($sample,clip) $clip] \
			-singlecell [get optionsa($sample,singlecell) $singlecell] \
			-sc_whitelist [get optionsa($sample,sc_whitelist) $singlecell_whitelist] \
			-sc_umisize [get optionsa($sample,sc_umisize) $singlecell_umisize] \
			-sc_filters [get optionsa($sample,sc_filters) $sc_filters] \
			-sc_celltypers [get optionsa($sample,sc_celltypers) $sc_celltypers] \
			-sc_expectedcells [get optionsa($sample,sc_expectedcells) $sc_expectedcells] \
			-cellmarkerfile [get optionsa($sample,cellmarkerfile) $cellmarkerfile] \
			-tissue [get optionsa($sample,tissue) $tissue] \
			-datatype [get optionsa($sample,datatype) $datatype] \
			-aliformat [get optionsa($sample,aliformat) $aliformat] \
			-aligners [get optionsa($sample,aligners) $aligners] \
			-ali_keepcomments [get optionsa($sample,ali_keepcomments) $ali_keepcomments] \
			-realign [get optionsa($sample,realign) $realign] \
			-varcallers [get optionsa($sample,varcallers) $varcallers] \
			-svcallers [get optionsa($sample,svcallers) $svcallers] \
			-methcallers [get optionsa($sample,methcallers) $methcallers] \
			-counters [get optionsa($sample,counters) $counters] \
			-reftranscripts [get optionsa($sample,reftranscripts) $reftranscripts] \
			-isocallers [get optionsa($sample,isocallers) $isocallers] \
			-hap_bam [get optionsa($sample,hap_bam) $hap_bam] \
			-dbdir [get optionsa($sample,dbdir) $dbdir] \
			-split [get optionsa($sample,split) $split] \
			-paired [get optionsa($sample,paired) $paired] \
			-maxfastqdistr [get optionsa($sample,maxfastqdistr) $maxfastqdistr] \
			-adapterfile [get optionsa($sample,adapterfile) $adapterfile] \
			-reports [get optionsa($sample,reports) $reports] \
			-samBQ [get optionsa($sample,samBQ) $samBQ] \
			-cleanup [get optionsa($sample,cleanup) $cleanup] \
			-removeduplicates [get optionsa($sample,removeduplicates) $removeduplicates] \
			-amplicons [get optionsa($sample,amplicons) $amplicons] \
			-threads [get optionsa($sample,threads) $threads] \
			-distrreg [get optionsa($sample,distrreg) $distrreg] \
			-keepsams [get optionsa($sample,keepsams) $keepsams] \
			-removeskew [get optionsa($sample,removeskew) $removeskew] \
			-dt [get optionsa($sample,dt) $dt] \
			-targetfile [get optionsa($sample,targetfile) $targetfile] \
			-minfastqreads [get optionsa($sample,minfastqreads) $minfastqreads] \
			-depth_histo_max [get optionsa($sample,depth_histo_max) $depth_histo_max] \
			$dir \
		]
		if {!$jobsample} {
			process_sample_job -todoVar todo {*}$sampleargs
		} else {
			# find deps and targets by running the process_sample_job with job_getinfo set to 1
			job_getinfo 1
			set verbose [logverbose]
			set ::deps {} ; set ::targets {}
			process_sample_job -todoVar todo {*}$sampleargs
			foreach {deps targets} [job_getinfo 0] break
			logverbose $verbose
			# run the actual job with deps and targets found
			job process_sample-$sample -deps $deps -targets $targets -vars {
				sampleargs
			} -code {
				cg process_sample -stack 1 -v 2 {*}$sampleargs >@ stdout 2>@ stderr
			}
		}
	}
	# run msamples (not supporting jobsample for these)
	foreach sample $msamples {
		putslog "Processing msample $sample"
		set dir $destdir/msamples/$sample
		process_sample_job \
			-preset [get optionsa($sample,preset) $preset] \
			-clip [get optionsa($sample,clip) $clip] \
			-singlecell [get optionsa($sample,singlecell) $singlecell] \
			-sc_whitelist [get optionsa($sample,sc_whitelist) $singlecell_whitelist] \
			-sc_umisize [get optionsa($sample,sc_umisize) $singlecell_umisize] \
			-sc_filters [get optionsa($sample,sc_filters) $sc_filters] \
			-sc_celltypers [get optionsa($sample,sc_celltypers) $sc_celltypers] \
			-sc_expectedcells [get optionsa($sample,sc_expectedcells) $sc_expectedcells] \
			-cellmarkerfile [get optionsa($sample,cellmarkerfile) $cellmarkerfile] \
			-tissue [get optionsa($sample,tissue) $tissue] \
			-datatype [get optionsa($sample,datatype) $datatype] \
			-aliformat [get optionsa($sample,aliformat) $aliformat] \
			-aligners [get optionsa($sample,aligners) $aligners] \
			-ali_keepcomments [get optionsa($sample,ali_keepcomments) $ali_keepcomments] \
			-realign [get optionsa($sample,realign) $realign] \
			-varcallers [get optionsa($sample,varcallers) $varcallers] \
			-svcallers [get optionsa($sample,svcallers) $svcallers] \
			-methcallers [get optionsa($sample,methcallers) $methcallers] \
			-counters [get optionsa($sample,counters) $counters] \
			-reftranscripts [get optionsa($sample,reftranscripts) $reftranscripts] \
			-isocallers [get optionsa($sample,isocallers) $isocallers] \
			-hap_bam [get optionsa($sample,hap_bam) $hap_bam] \
			-dbdir [get optionsa($sample,dbdir) $dbdir] \
			-split [get optionsa($sample,split) $split] \
			-paired [get optionsa($sample,paired) $paired] \
			-maxfastqdistr [get optionsa($sample,maxfastqdistr) $maxfastqdistr] \
			-adapterfile [get optionsa($sample,adapterfile) $adapterfile] \
			-reports [get optionsa($sample,reports) $reports] \
			-samBQ [get optionsa($sample,samBQ) $samBQ] \
			-cleanup [get optionsa($sample,cleanup) $cleanup] \
			-removeduplicates [get optionsa($sample,removeduplicates) $removeduplicates] \
			-amplicons [get optionsa($sample,amplicons) $amplicons] \
			-threads [get optionsa($sample,threads) $threads] \
			-distrreg [get optionsa($sample,distrreg) $distrreg] \
			-keepsams [get optionsa($sample,keepsams) $keepsams] \
			-removeskew [get optionsa($sample,removeskew) $removeskew] \
			-dt [get optionsa($sample,dt) $dt] \
			-targetfile [get optionsa($sample,targetfile) $targetfile] \
			-minfastqreads [get optionsa($sample,minfastqreads) $minfastqreads] \
			-depth_histo_max [get optionsa($sample,depth_histo_max) $depth_histo_max] \
			$dir
	}
	set_job_logdir $destdir/log_jobs
	set todo(var) [list_remdup $todo(var)]
	set todo(sv) [list_remdup $todo(sv)]
	set todo(meth) [list_remdup $todo(meth)]
	set todo(reports) [list_remdup $todo(reports)]
	if {![llength $reports]} {set todo(reports) {}}
	process_multicompar_job \
		-experiment $experiment \
		-skipincomplete 1 -targetvarsfile $targetvarsfile \
		-varfiles $todo(var) -svfiles $todo(sv) -methfiles $todo(meth) \
		-counters $counters \
		-isocallers $isocallers \
		-iso_match $iso_match \
		-sc_celltypers $sc_celltypers \
		-threads $threads \
		-distrreg $distrreg \
		-keepfields $keepfields \
		-split $split -dbfiles $dbfiles \
		-cleanup $cleanup \
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
