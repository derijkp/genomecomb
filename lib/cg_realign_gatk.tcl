proc realign_gatk_job {args} {
	upvar job_logdir job_logdir
	set skips {}
	set regionfile {}
	set threads 2
	cg_options realign_gatk args {
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
		realign around indels using gatk
	}
	set bamfile [file_absolute $bamfile]
	set resultbamfile [file_absolute $resultbamfile]
	set refseq [refseq $refseq]
	if {[file isdir $refseq]} {
		set refseq [lindex [glob $refseq/genome_*.ifas] 0]
	}
	if {![info exists job_logdir]} {
		job_logdir $resultbamfile.log_jobs
	}
	if {$regionfile eq ""} {
		set cov3reg [bam2reg_job -mincoverage 3 $bamfile]
		set regionfile $cov3reg
	}
	set gatkrefseq [gatk_refseq_job $refseq]
	set dict [file root $gatkrefseq].dict
	job realign_gatk-[file tail $resultbamfile] -mem [job_mempercore 10G [expr {1+$threads}]] -cores [expr {1+$threads}] \
	-deps {$bamfile ($$bamfile.bai) $dict $gatkrefseq $refseq $regionfile} \
	-targets {$resultbamfile $resultbamfile.bai} {*}$skips \
	-vars {gatkrefseq refseq gatk bamfile regionfile threads} -code {
		putslog "making $target"
		analysisinfo_write $dep $target realign gatk realign_version [version gatk3]
		if {![file exists $bamfile.bai]} {exec samtools index $bamfile}
		set bedfile [tempbed $regionfile $refseq]
		lappend realignopts -L $bedfile
		gatk3exec {-XX:ParallelGCThreads=1 -Xms512m -Xmx8g} RealignerTargetCreator -R $gatkrefseq -I $dep -o $target.intervals {*}$realignopts 2>@ stdout >@ stdout
		if {[loc_compare [version gatk3] 2.7] >= 0} {
			set extra {--filter_bases_not_stored}
		} else {
			set extra {}
		}
		lappend extra --filter_mismatching_base_and_quals
		gatk3exec {-XX:ParallelGCThreads=1 -Xms512m -Xmx8g} IndelRealigner -R $gatkrefseq \
			-targetIntervals $target.intervals -I $dep \
			-o $target.temp {*}$extra 2>@ stdout >@ stdout
		catch {file rename -force $target.temp.bai $target.bai}
		catch {file delete $target.intervals}
		file rename -force $target.temp $target
	}
}

proc cg_realign_gatk {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [realign_gatk_job {*}$args]
	job_wait
	return $result
}
