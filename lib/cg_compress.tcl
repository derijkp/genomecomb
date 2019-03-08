proc cg_compress_job args {
	set cmdline [list cg compress {*}$args]
	set keep {}
	set compressionlevel [defcompressionlevel]
	set blocksize 5
	set index 0
	set threads 1
	set outputfile {}
	set method lz4
	cg_options compress args {
		-o - -outputfile {
			set outputfile $value
			if {$keep eq ""} {set keep 1}
		}
		-k - -keep {
			set keep $value
		}
		-i - -index {
			set index $value
		}
		-t - -threads {
			set threads $value
		}
		-c - -compressionlevel {
			set compressionlevel $value
		}
		-b - -blocksize {
			set blocksize $value
		}
		-m - -method {
			set method $value
		}
	}
	if {$keep eq ""} {set keep 0}
	if {$outputfile ne "" && [llength $args] > 1} {
		error "option -o can only be used for compressing one file"
	}
	if {![llength $args]} {
		set args [list -]
		if {$outputfile eq ""} {set outputfile -}
	}
	if {[job_distribute] ne "0"} {
		upvar job_logdir job_logdir
		if {$args eq "-" || $outputfile eq "-"} {
			error "cannot run compression job distributed on stdin/stdout"
		}
		set logfile [job_logfile compress_$method [pwd] $cmdline {*}[versions $method]]
	}
	foreach file $args {
		if {$outputfile eq ""} {
			set target [file root $file].$method
		} else {
			set target $outputfile
		}
		if {[job_distribute] eq "0"} {
			compress_$method $file $target $index $keep $threads $compressionlevel $blocksize
		} else {
			job $method-$file -deps {$file} -targets {$target} -vars {file keep index threads compressionlevel blocksize} -code {
				compress_$method $file $target $index $keep $threads $compressionlevel $blocksize
			}
		}
	}
}

proc cg_compress args {
	set args [job_init {*}$args]
	cg_compress_job {*}$args
	job_wait
}
