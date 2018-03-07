proc realign_abra_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	set regionfile {}
	set threads 2
	cg_options map_minimap2 args {
		-skips {
			set skips $value
		}
		-regionfile {
			set regionfile $value
		}
		-threads - -t {
			set threads $value
		}
	} {bamfile resultbamfile refseq} 3 3 {
		realign around indels using abra
	}
	set bamfile [file_absolute $bamfile]
	set resultbamfile [file_absolute $resultbamfile]
	set refseq [file_absolute $refseq]
	if {[file isdir $refseq]} {
		set refseq [lindex [glob $refseq/genome_*.ifas] 0]
	}
	if {![info exists job_logdir]} {
		job_logdir $resultbamfile.log_jobs
	}
	set abra [findjar abra]
	job realign_abra-[file tail $resultbamfile] -cores $threads \
	-deps {$bamfile $refseq ($regionfile)} \
	-targets {$resultbamfile $resultbamfile.bai} {*}$skips \
	-vars {abra bamfile refseq regionfile threads} -code {
		if {$regionfile eq "" || ![file exists $regionfile]} {
			putslog "making regionfile"
			set regionfile [tempfile]
			cg bam2reg $bamfile 3 $regionfile
		}
		putslog "making $target"
		if {![file exists $bamfile.bai]} {exec samtools index $bamfile}
		if {[file extension $regionfile] ne ".bed"} {
			set tempfile [tempfile]
			cg tsv2bed $regionfile $tempfile
			set regionfile $tempfile
		}
		if {[catch {
			exec java -Xmx4G -XX:ParallelGCThreads=1 -jar $abra --in $bamfile --out $target.temp.bam --ref $refseq --targets $regionfile --threads $threads --working [scratchdir] > $target.log 2>@ stdout 1>@ stdout
		} msg]} {
			error $msg
		}
		exec samtools index $target.temp.bam
		file rename -force $target.temp.bam.bai $target.bai
		file rename -force $target.temp.bam $target
	}
}

proc cg_realign_abra {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [realign_abra_job {*}$args]
	job_wait
	return $result
}
