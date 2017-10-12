proc basecaller_albacore_job {args} {
	upvar job_logdir job_logdir
	set keepargs $args
	set reads 50000
	set config {}
	set flowcell {}
	set kit {}
	set opts {}
	set threads 1
	set barcoding 0
	cg_options var_sam args {
		-r - -reads {
			set reads $value
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
	} {sourcedir resultdir}
	set sourcedir [file_absolute $sourcedir]
	set resultdir [file_absolute $resultdir]
	if {$config eq ""} {
		if {$kit eq "" || $flowcell eq ""} {
			error "Either -config or -flowcell and -kit must be specified"
		}
		lappend opts -kit $kit -flowcell $flowcell
	} else {
		lappend opts --config $config
	}
	job_logfile $resultdir/basecall_albacore_[file tail $resultdir] $resultdir \
		[list cg basecaller_albacore {*}$keepargs] \
		{*}[versions albacore]
	# find fast5 files
	set files [lsort -dict [dirglob $sourcedir *.fast5]]
	set len [llength $files]
	file mkdir $resultdir
	file_write $resultdir/info.txt "source\t$sourcedir\nsize\t$len\n"
	for {set from 0} {$from < $len} {incr from $reads} {
		set to [expr {$from+$reads-1}]
		if {$to > $len} {set to $len}
		set todo [lrange $files $from $to]
		set target1 $resultdir/pass_[file tail $resultdir]-$from-$to.fastq.gz
		set target2 $resultdir/fail_[file tail $resultdir]-$from-$to.fastq.gz
		job albacore-[file tail $sourcedir]-$from-$to -cores $threads \
		-deps $sourcedir \
		-targets {$target1 $target2} \
		-vars {sourcedir from to todo opts threads} \
		-code {
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

proc cg_basecaller_albacore {args} {
	set args [job_init {*}$args]
	basecaller_albacore_job {*}$args
	job_wait
}
