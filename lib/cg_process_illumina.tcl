proc cg_bcl2fastq {rundir outdir {rtr 6} {dtr 6} {ptr 6} {wtr 6} } {
	#-r, --loading-threads Number of threads used for loading BCL data.
	#-d, --demultiplexing-threads Number of threads used for demultiplexing.
	#-p, --processing-threads  Number of threads used for processing demultiplexed data.
	#-w, --writing-threads Number of threads used for writing FASTQ data. This must not be higher than number of samples.
	exec [bcl2fastq] --create-fastq-for-index-reads -r $rtr -d $dtr -p $ptr -w $wtr --runfolder-dir $rundir --output-dir $outdir 2>@ stderr >@ stdout
}

proc cg_process_conv_illnextseq {illsrc destdir} {
	set illsrc [file_absolute $illsrc]
	set destdir [file_absolute $destdir]
	file mkdir $destdir
	set keeppwd [pwd]
	cd $destdir
	# copy files from illumina dir
	set files [glob $illsrc/*_R*.fastq*]
	foreach file $files {
		set sample [file tail $file]
		regsub {_[^_]+_[^_]+_[^_]+_[^_]+\.fastq.*} $sample {} sample
		regsub -all -- - $sample _ sample 
		if {$sample != "Undetermined"} {
			file mkdir $destdir/$sample/ori
			file mkdir $destdir/$sample/ori/fastq
			file mkdir $destdir/$sample
			file mkdir $destdir/$sample/fastq
			hardlink $file $destdir/$sample/ori/fastq
			cplinked $destdir/$sample/ori/fastq $destdir/$sample/fastq
		}
	}
	cd $keeppwd
}

proc process_illumina {args} {
	set dbdir {}
	set dbfiles {}
	set targetfile {}
	set realign 1
	set cleanup 1
	set cleanupfiles {}
	set paired 1
	set adapterfile {}
	set conv_nextseq 0
	set reports all
	cg_options process_illumina args {
		-realign {
			set realign $value
		}
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
		-targetfile {
			set targetfile [file_absolute $value]
			if {$value ne "" && ![jobfileexists $targetfile]} {error "target file $targetfile does not exists"}
		}
		-paired {
			set paired $value
		}
		-adapterfile {
			if {$value ne "" && ![file exists $value]} {error "adapterfile $value does not exists"}
			set adapterfile [file_absolute $value]
		}
		-conv_nextseq {
			set conv_nextseq $value
		}
		-reports {
			set reports $value
		}
		-c - -cleanup {
			set cleanup $value
		}
		-m - -maxopenfiles {
			set maxopenfiles $value
			set ::maxopenfiles [expr {$maxopenfiles - 4}]
		}
	} {destdir dbdir} 1 2
	set destdir [file_absolute $destdir]
	set dbdir [file_absolute $dbdir]
	set adapterfile [adapterfile $adapterfile]
	set experimentname [file tail $destdir]
	# check projectinfo
	projectinfo $destdir dbdir {split 1}
	set dbdir [dbdir $dbdir]
	# logfile
	set cmdline [list cg process_illumina]
	foreach option {
		realign dbdir dbfiles paired adapterfile conv_nextseq reports cleanup maxopenfiles
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline $destdir $dbdir
	job_logfile $destdir/process_illumina_[file tail $destdir] $destdir $cmdline \
		{*}[versions dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk biobambam picard java gnusort8 lz4 os]
	# start
	##in case of nextseq500 data - generate fastqs & distribute data
	if {$conv_nextseq} {
		set rundir [glob $destdir/*NS500*]
		cg_bcl2fastq $rundir fastq 4 4 4 4
		cg_process_conv_illnextseq fastq $destdir
	}
	set refseq [glob $dbdir/genome_*.ifas]
	set resultbamprefix rds
	set samples {}
	set experiment [file tail $destdir]
	if {[file exists $destdir/samples]} {
		set sampledir $destdir/samples
	} else {
		set sampledir $destdir
	}
	foreach dir [dirglob $sampledir */fastq] {
		lappend samples [file dir $dir]
	}
	set samples [ssort -natural $samples]
	set keeppwd [pwd]
	cd $destdir
	job_logdir $destdir/log_jobs
	set todo {}
	set reportstodo {}
	foreach sample $samples {
		puts $sample
		set dir $sampledir/$sample
		catch {file mkdir $dir}
		puts $dir
		cd $dir
		job_logdir $dir/log_jobs
		# find fastq files
		set fastqfiles [ssort -natural [jobglob fastq/*.fastq.gz fastq/*.fastq fastq/*.fq.gz fastq/*.fq]]
		if {![llength $fastqfiles]} {
			# if there are no fastqfiles, check if there are bam files in ori to extract fastq from
			set files [ssort -natural [jobglob ori/*.bam]]
			foreach file $files {
				set base fastq/[file tail [file root $file]]
				set target $base-R1.fastq.gz
				set target2 $base-R2.fastq.gz
				job bam2fastq-[file tail $file] -deps {$file} \
				-targets {$target $target2} -code {
					cg bam2fastq $dep $target.temp.gz $target2.temp.gz
					file rename -force $target.temp.gz $target
					file rename -force $target2.temp.gz $target2
				}
			}
			set fastqfiles [ssort -natural [jobglob fastq/*.fastq.gz fastq/*.fastq fastq/*.fq.gz fastq/*.fq]]
		}
		set resultbamfile map-${resultbamprefix}bwa-$sample.bam
		if {[llength $fastqfiles]} {
			# do not do any of preliminaries if end product is already there
			set bamfile map-sbwa-$sample.bam
			# quality and adapter clipping
			set files [fastq_clipadapters_job $fastqfiles \
				-adapterfile $adapterfile -paired $paired \
				-skips [list -skip [list $bamfile $bamfile.analysisinfo] -skip [list $resultbamfile $resultbamfile.analysisinfo]]]
			lappend cleanupfiles {*}$files [file dir [lindex $files 0]]
			# map using bwa
			map_bwa_job -paired $paired -skips [list -skip [list $resultbamfile $resultbamfile.analysisinfo]] $bamfile $refseq $sample {*}$files
		}
		# extract regions with coverage >= 5 (for cleaning)
		set cov5reg [bam2reg_job -mincoverage 5 -skip [list $resultbamfile $resultbamfile.analysisinfo] map-sbwa-$sample.bam]
		# clean bamfile (mark duplicates, realign)
		set cleanedbam [bam_clean_job \
			-sort 0 -removeduplicates 1 -realign $realign -regionfile $cov5reg -cleanup $cleanup \
			 map-sbwa-$sample.bam $refseq $sample]
		# make 5x coverage regfile from cleanedbam
		set cov5reg [bam2reg_job -mincoverage 5 $cleanedbam]
		# make 20x coverage regfile
		bam2reg_job -mincoverage 20 -compress 1 $cleanedbam
		# gatk variant calling on map-rdsbwa
		var_gatk_job -regionfile $cov5reg -split $split $cleanedbam $refseq
		lappend todo gatk-rdsbwa-$sample
		# samtools variant calling on map-rdsbwa
		var_sam_job -regionfile $cov5reg -split $split $cleanedbam $refseq
		lappend todo sam-rdsbwa-$sample
		if {$cleanup} {
			# clean up no longer needed intermediate files
			cleanup_job cleanupsample-$sample $cleanupfiles [list $resultbamfile var-gatk-rdsbwa-$sample.tsv var-sam-rdsbwa-$sample.tsv]
		}
		# convert existing vcfs
		set files [ssort -natural [jobglob var-*.vcf]]
		foreach file $files {
			set target [file root [gzroot $file]].tsv
			if {![file exists $target]} {
				job vcf2tsv-$file -deps {$file} -targets {$target} -vars split -code {
					cg vcf2tsv -split $split $dep $target.temp
					file rename -force $target.temp $target
				}
				lappend todo [string range $target 4 end-4]
			}
		}
		# add existing var files to todo
		# This may add duplicates (already existing var files from previous gatk,sam), but will be cleaned later by he list_remdup
		set files [ssort -natural [jobglob var-*.tsv]]
		foreach file $files {
			set target [file root [gzroot $file]].tsv
			lappend todo [string range $target 4 end-4]
		}
		# calculate reports
		if {[llength $reports]} {
			process_reports_job $sampledir/$sample $dbdir $reports
			lappend reportstodo $sampledir/$sample/reports
		}
	}
	job_logdir $destdir/log_jobs
	cd $destdir
	set todo [list_remdup $todo]
	process_multicompar_job -experiment $experiment -skipincomplete 1 -split $split -dbfiles $dbfiles $destdir $dbdir $todo
	if {[llength $reports]} {
		proces_reportscombine_job $destdir/reports {*}$reportstodo
		if {[jobfileexists $destdir/reports/report_hsmetrics-${experimentname}.tsv]} {
			mklink $destdir/reports/report_hsmetrics-${experimentname}.tsv $destdir/${experimentname}_hsmetrics_report.tsv
		}
	}
	cd $keeppwd
}

proc cg_process_illumina {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_illumina
	}
	process_illumina {*}$args
	job_wait
}
