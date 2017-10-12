proc minimap2refseq_job {refseq preset} {
	upvar job_logdir job_logdir
	set minimap2refseq $refseq.minimap2.$preset
	if {[file exists $minimap2refseq]} {return $minimap2refseq}
	set tail [file tail $refseq]
	if {[jobtargetexists [list $minimap2refseq] $refseq]} return
	job minimap2_2refseq-[file tail $refseq] -deps {$refseq} -targets {$minimap2refseq} -vars {preset} -code {
		if {[catch {exec -ignorestderr minimap2 -x $preset -d $target $dep} e]} {
			error $e
		}
	}
	return $minimap2refseq
}

proc map_minimap2_job {args} {
	upvar job_logdir job_logdir
	set preset map-ont
	set readgroupdata {}
	set skips {}
	set threads 2
	cg_options map_minimap2 args {
		-x - -preset - -p {
			set preset $value
		}
		-readgroupdata {
			set readgroupdata $value
		}
		-skips {
			set skips $value
		}
		-threads - -t {
			set threads $value
		}
	} {result refseq sample fastqfile1} 4
	set files [list $fastqfile1 {*}$args]
	if {![info exists job_logdir]} {
		job_logdir [file dir $result]/log_jobs
	}
	array set a [list PL illumina LB solexa-123 PU $sample SM $sample]
	if {$readgroupdata ne ""} {
		array set a $readgroupdata
	}
	set readgroupdata [array get a]
	set minimap2refseq [minimap2refseq_job $refseq $preset]
	set resultbase [file root $result]
	set samfiles {}
	set num 1
	foreach file $files {
		set name [file root [file tail $file]]
		set target $resultbase-$name.sam
		lappend samfiles $target
		job minimap2-$sample-$name -mem 5G -cores $threads \
		-deps [list $minimap2refseq $file] -targets {$target} \
		-vars {threads preset readgroupdata sample} \
		-skip [list $resultbase.bam] {*}$skips -code {
			puts "making $target"
			foreach {minimap2refseq fastq} $deps break
			set rg {}
			foreach {key value} $readgroupdata {
				lappend rg "$key:$value"
			}
			exec minimap2 -a -x $preset -t $threads -R @RG\\tID:$sample\\t[join $rg \\t] \
				$minimap2refseq $fastq > $target.temp 2>@ stderr
			file rename -force $target.temp $target
		}
	}
	job minimap2_2bam-$sample -deps $samfiles -rmtargets $samfiles -targets {$result $result.bai} {*}$skips -vars {resultbase} -code {
		puts "making $target"
		if {[catch {
			exec samcat {*}$deps | bamsort SO=coordinate tmpfile=[scratchfile] index=1 indexfilename=$target.bai inputformat=sam > $target.temp 2>@ stderr
		}]} {
			error $msg
		}
		file rename -force $target.temp $target
		file delete {*}$deps
	}
}

proc cg_map_minimap2 {args} {
	set args [job_init {*}$args]
	unset job_logdir
	map_minimap2_job {*}$args
	job_wait
}
