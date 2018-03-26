proc minimap2refseq_job {refseq preset} {
	upvar job_logdir job_logdir
	set minimap2refseq $refseq.minimap2.$preset
	if {[file exists $minimap2refseq]} {return $minimap2refseq}
	set tail [file tail $refseq]
	if {[jobtargetexists [list $minimap2refseq] $refseq]} return
	job minimap2_2refseq-[file tail $refseq] -deps {$refseq} -targets {$minimap2refseq} -vars {preset} -code {
		if {[catch {exec -ignorestderr minimap2 -x $preset -d $target.temp $dep} e]} {
			error $e
		}
		file rename $target.temp $target
	}
	return $minimap2refseq
}

proc map_minimap2_job {args} {
	upvar job_logdir job_logdir
	set paired 1
	set keepargs $args
	set preset {}
	set readgroupdata {}
	set skips {}
	set keepsams 0
	set threads 2
	cg_options map_minimap2 args {
		-paired - -p {
			set paired $value
		}
		-x - -preset - -p {
			set preset $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-keepsams {
			set keepsams [true $value]
		}
		-skips {
			set skips $value
		}
		-threads - -t {
			set threads $value
		}
	} {result refseq sample fastqfile1} 4 ... {
		align reads in fastq files to a reference genome using minimap2
	}
	if {$paired && $preset eq ""} {set preset sr} else {set preset map-ont}
	set files [list $fastqfile1 {*}$args]
	set result [file_absolute $result]
	set refseq [file_absolute $refseq]
	set resultdir [file dir $result]
	file mkdir $result.index
	if {![info exists job_logdir]} {
		job_logdir $result.index/log_jobs
	}
	job_logfile $resultdir/map_minmap2_[file tail $result] $resultdir \
		[list cg map_minmap2_ {*}$keepargs] \
		{*}[versions minimap2]
	#
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set readgroupdata [array get a]
	set minimap2refseq [minimap2refseq_job $refseq $preset]
	set resultbase [file root $result]
	set samfiles {}
	set num 1
	if {!$paired} {
		foreach file $files {
			set name [file root [file tail $file]]
			set target $result.index/$name.sam
			lappend samfiles $target
			set analysisinfo [gzroot $target].analysisinfo
			job minimap2-$sample-$name -mem 5G -cores $threads -deps {
				$minimap2refseq $file
			} -targets {
				$target $analysisinfo
			} -vars {
				threads preset readgroupdata sample
			} -skip [list $resultbase.bam] {*}$skips -code {
				puts "making $target"
				foreach {minimap2refseq fastq} $deps break
				analysisinfo_write $fastq $target aligner minimap2 aligner_version [version minimap2] reference [file2refname $minimap2refseq] aligner_paired 0
				set rg {}
				foreach {key value} $readgroupdata {
					lappend rg "$key:$value"
				}
				exec minimap2 -a -x $preset -t $threads -R @RG\\tID:$sample\\t[join $rg \\t] \
					$minimap2refseq $fastq > $target.temp 2>@ stderr
				file rename -force $target.temp $target
			}
		}
	} else {
		foreach {file1 file2} $files {
			set name [file root [file tail $file1]]
			set target $resultbase-$name.sam
			lappend samfiles $target
			set analysisinfo [gzroot $target].analysisinfo
			lappend asamfiles $analysisinfo
			job bwa-$sample-$name -mem 5G -cores $threads -deps {
				$minimap2refseq $file1 $file2
			} -targets {
				$target $analysisinfo
			} -vars {
				threads preset readgroupdata sample
			} -skip [list $resultbase.bam] {*}$skips -code {
				puts "making $target"
				foreach {minimap2refseq fastq1 fastq2} $deps break
				analysisinfo_write $fastq1 $target aligner minimap2 aligner_version [version minimap2] reference [file2refname $minimap2refseq] aligner_paired 1
				set rg {}
				foreach {key value} $readgroupdata {
					lappend rg "$key:$value"
				}
				exec minimap2 -a -x $preset -t $threads -R @RG\\tID:$sample\\t[join $rg \\t] \
					$minimap2refseq $fastq1 $fastq2 > $target.temp 2>@ stderr
				file rename -force $target.temp $target
			}
		}
	}
	if {$keepsams} {set rmsamfiles {}} else {set rmsamfiles $samfiles}
	job minimap2_2bam-$sample -deps $samfiles -rmtargets $rmsamfiles \
	-targets {$result $result.bai $result.analysisinfo} {*}$skips -vars {resultbase} -code {
		puts "making $target"
		analysisinfo_write $dep $target
		if {[catch {
			exec samcat {*}$deps | bamsort SO=coordinate tmpfile=[scratchfile] index=1 indexfilename=$target.bai inputformat=sam > $target.temp 2>@ stderr
		}]} {
			error $msg
		}
		file rename -force $target.temp $target
		foreach dep $deps {
			file delete $dep [gzroot $dep].analysisinfo
		}
	}
}

proc cg_map_minimap2 {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [map_minimap2_job {*}$args]
	job_wait
	return $result
}
