proc bwarefseq_job {refseq} {
	upvar job_logdir job_logdir
	set bwarefseq $refseq.bwa/[file tail $refseq]
	set bwarefseqfa [file root $bwarefseq].fa
	if {[file exists $bwarefseqfa]} {return $bwarefseqfa}
	set tail [file tail $refseq]
	set targets [list $refseq.bwa $refseq.bwa/$tail $refseq.bwa/$tail.fai]
	foreach ext {amb ann bwt pac sa} {lappend targets $refseq.bwa/$tail.$ext}
	if {[jobtargetexists $targets $refseq]} return
	job bwa2refseq-[file tail $refseq] -deps {$refseq} -targets $targets -code {
		file delete -force $target.temp
		file mkdir $target.temp
		mklink $dep $target.temp/[file tail $dep]
		if {![file exists $dep.fai]} {
			exec samtools faidx $dep
		}
		mklink $dep.fai $target.temp/[file tail $dep].fai
		exec bwa index $target.temp/[file tail $dep] 2>@ stderr
		file delete -force $target
		file rename $target.temp $target
	}
	# with .fa as extension
	set deps [jobglob $refseq.bwa/*.ifas*]
	set targets {}
	foreach file $deps {
		if {[file extension $file] eq ".ifas"} {
			lappend targets [file root $file].fa
		} else {
			if {![regsub {\.ifas\.([^.]+)$} $file {.fa.\1} nfile]} {error "could not convert name of $file"}
			lappend targets $nfile
		}
	}
	job bwa2refseq_fa-[file tail $refseq] -deps $deps -targets $targets -code {
		foreach file $deps nfile $targets {
			mklink $file $nfile
		}
	}
	return $bwarefseqfa
}

proc map_bwa_job {args} {
	upvar job_logdir job_logdir
	set paired 1
	set readgroupdata {}
	set pre {}
	set skips {}
	set threads 2
	set keepsams 0
	cg_options map_bwa args {
		-paired - -p {
			set paired $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-pre {
			set pre $value
		}
		-skips {
			set skips $value
		}
		-keepsams {
			set keepsams $value
		}
		-threads - -t {
			set threads $value
		}
	} {result refseq sample fastqfile1} 4 ... {
		align reads in fastq files to a reference genome using bwa-mem
	}
	set files [list $fastqfile1 {*}$args]
	if {![info exists job_logdir]} {
		job_logdir [file dir $result]/log_jobs
	}
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set readgroupdata [array get a]
	set bwarefseq [bwarefseq_job $refseq]
	set resultbase [file root $result]
	set samfiles {}
	set asamfiles {}
	set num 1
	if {!$paired} {
		foreach file $files {
			set name [file root [file tail $file]]
			set target $resultbase-$name.sam
			lappend samfiles $target
			set analysisinfo [gzroot $target].analysisinfo
			lappend asamfiles $analysisinfo
			job bwa-$sample-$name -mem 5G -cores $threads \
			-deps [list $bwarefseq $file] \
			-targets {$target $analysisinfo} \
			-vars {readgroupdata sample paired threads} \
			-skip [list $resultbase.bam] {*}$skips -code {
				puts "making $target"
				foreach {bwarefseq fastq} $deps break
				analysisinfo_write $fastq $target aligner bwa aligner_version [version bwa] reference [file2refname $bwarefseq] aligner_paired 0
				set rg {}
				foreach {key value} $readgroupdata {
					lappend rg "$key:$value"
				}
				exec bwa mem -t $threads -R @RG\\tID:$sample\\t[join $rg \\t] $bwarefseq $fastq > $target.temp 2>@ stderr
				file rename -force $target.temp $target
			}
		}
	} else {
		if {[expr {[llength $files]%2}]} {
			error "bwa needs even number of files for paired analysis"
		}
		foreach {file1 file2} $files {
			set name [file root [file tail $file1]]
			set target $resultbase-$name.sam
			lappend samfiles $target
			set analysisinfo [gzroot $target].analysisinfo
			lappend asamfiles $analysisinfo
			job bwa-$sample-$name -mem 5G -cores $threads \
			-deps [list $bwarefseq $file1 $file2] \
			-targets {$target $analysisinfo} \
			-vars {readgroupdata sample paired threads} \
			-skip [list $resultbase.bam] {*}$skips -code {
				puts "making $target"
				foreach {bwarefseq fastq1 fastq2} $deps break
				analysisinfo_write $fastq1 $target aligner bwa aligner_version [version bwa] reference [file2refname $bwarefseq] aligner_paired 1
				set rg {}
				foreach {key value} $readgroupdata {
					lappend rg "$key:$value"
				}
				exec bwa mem -t $threads -M -R @RG\\tID:$sample\\t[join $rg \\t] $bwarefseq $fastq1 $fastq2 > $target.temp 2>@ stderr
				file rename -force $target.temp $target
				file delete bwa1.fastq bwa2.fastq
			}
		}
	}
	sam_catmerge_job -skips $skips -name bwa2bam-$sample -deletesams [string is false $keepsams] -threads $threads $result {*}$samfiles
}

proc cg_map_bwa {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set return [map_bwa_job {*}$args]
	job_wait
	return $return
}
