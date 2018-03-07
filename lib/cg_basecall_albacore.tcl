proc basecaller_albacore_job {args} {
	upvar job_logdir job_logdir
	set keepargs $args
	set numreads 4000
	set subdirs 1
	set config {}
	set flowcell {}
	set kit {}
	set opts {}
	set threads 1
	set barcoding 0
	cg_options var_sam args {
		-n - -numreads {
			set numreads $value
		}
		-s - -subdirs {
			set subdirs $value
		}
		-c - -config {
			set config $value
		}
		-f - -flowcell {
			set flowcell $value
		}
		-k - -kit {
			set kit $value
		}
		-filter {
			if {![true $value]} {lappend opts --disable_filtering}
		}
		-b - -barcoding {
			set barcoding [true $value]
			if {$barcoding} {
				lappend opts --barcoding
			}
		}
		-t - -threads {
			set threads $value
		}
		default {
			lappend opts $key $value
		}
	} {resultdir sourcedir} 2 2 {
		basecall nanopore reads using albacore
	}
	set resultdir [file_absolute $resultdir]
	set sourcedirs [list [file_absolute $sourcedir]]
	foreach dir $args {
		lappend sourcedirs [file_absolute $dir]
	}
	if {$config eq ""} {
		if {$kit eq "" || $flowcell eq ""} {
			error "Either -config or -flowcell and -kit must be specified"
		}
		lappend opts -kit $kit -flowcell $flowcell
	} else {
		lappend opts --config $config
	}
	if {![info exists job_logdir]} {
		job_logdir $resultdir/log_jobs
	}
	set name [file tail $resultdir]
	if {$name eq "fastq"} {
		set name [file tail [file dir $resultdir]]
	}
	job_logfile $resultdir/basecall_albacore_$name $resultdir \
		[list cg basecaller_albacore {*}$keepargs] \
		{*}[versions albacore]

	if {$subdirs} {
		foreach sourcedir $sourcedirs {
			set sourcename [file tail $sourcedir]
			# find subdirs
			set dirs [lsort -dict [dirglob $sourcedir fast5/*]]
			file mkdir $resultdir
			file_write $resultdir/info.txt "source\t$sourcedir\nsubdirs\t[llength $dirs]\n"
			foreach dir $dirs {
				set files [dirglob $sourcedir/$dir *.fast5]
				set len [llength $files]
				if {$subdirs == 1 && $len < $numreads} {
					puts "skipping $dir: not full yet"
					continue
				}
				set name ${sourcename}_[file tail $dir]
				set target1 $resultdir/pass_$name.fastq.gz
				set target2 $resultdir/fail_$name.fastq.gz
				set target3 $resultdir/sequencing_summary_$name.txt
				
				job albacore-[file tail $sourcedir]-$name -cores $threads \
				-deps {$sourcedir/$dir} \
				-targets {$target1 $target2 $target3} \
				-vars {sourcedir from to todo opts threads} \
				-code {
					set tempdir [scratchdir]
					exec read_fast5_basecaller.py {*}$opts \
						--files_per_batch_folder 0 \
						--reads_per_fastq_batch 0 \
						-t $threads -o fastq \
						-i $dep \
						-s $tempdir/tempfastq >@ stdout 2>@ stderr
					exec cat {*}[glob $tempdir/tempfastq/workspace/pass/*.fastq] | gzip > $target1
					exec cat {*}[glob $tempdir/tempfastq/workspace/fail/*.fastq] | gzip > $target2
					file copy $tempdir/tempfastq/sequencing_summary.txt $target3
					file delete -force $tempdir/fast5 $tempdir/tempfastq
				}
			}
		}
	# old method of selecting 1 directory and finding all the *.fast5
	} else {
		# find fast5 files
		set files [lsort -dict [dirglob $sourcedir *.fast5]]
		set len [llength $files]
		if {$len == 0} {error "no fast5 files found (use -subdirs?)"}
		file mkdir $resultdir
		file_write $resultdir/info.txt "source\t$sourcedir\nsize\t$len\n"
		for {set from 0} {$from < $len} {incr from $numreads} {
			set to [expr {$from+$numreads-1}]
			if {$to > $len} {set to $len}
			set todo [lrange $files $from $to]
			set target1 $resultdir/pass_[file tail $resultdir]-$from-$to.fastq.gz
			set target2 $resultdir/fail_[file tail $resultdir]-$from-$to.fastq.gz
			job albacore-[file tail $sourcedir]-$from-$to -cores $threads \
			-deps $sourcedir \
			-targets {$target1 $target2} \
			-vars {sourcedir from to todo opts threads} \
			-code {
				albacoreinit
				set tempdir [scratchdir]
				file mkdir $tempdir/fast5
				foreach file $todo {
					mklink $sourcedir/$file $tempdir/fast5/$file 1
				}
				exec read_fast5_basecaller.py {*}$opts \
					--files_per_batch_folder 0 \
					--reads_per_fastq_batch 0 \
					-t $threads -o fastq \
					-i $tempdir/fast5 \
					-s $tempdir/tempfastq >@ stdout 2>@ stderr
				exec cat {*}[glob $tempdir/tempfastq/workspace/pass/*.fastq] | gzip > $target1
				exec cat {*}[glob $tempdir/tempfastq/workspace/fail/*.fastq] | gzip > $target2
				file delete -force $tempdir/fast5 $tempdir/tempfastq
			}
		}
	}

}

proc cg_basecaller_albacore {args} {
	set args [job_init {*}$args]
	set result [basecaller_albacore_job {*}$args]
	job_wait
	return $result
}
