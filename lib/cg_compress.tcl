proc cg_compress_job args {
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg compress {*}$args]
	set keep {}
	set compressionlevel {}
	set blocksize {}
	set index 0
	set threads {}
	set outputfile {}
	set method zst
	cg_options compress args {
		-o - -outputfile {
			set outputfile $value
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
			if {$value eq "gzip"} {set value gz}
			set method $value
		}
	}
	if {$outputfile ne "" && $keep eq ""} {set keep 1}
	if {$keep eq ""} {set keep 0}
	if {$outputfile ne "" && [llength $args] > 1} {
		error "option -o can only be used for compressing one file"
	}
	if {![llength $args]} {
		set args [list -]
		if {$outputfile eq ""} {set outputfile -}
	}
	if {[job_distribute] ne "0"} {
		if {$args eq "-" || $outputfile eq "-"} {
			error "cannot run compression job distributed on stdin/stdout"
		}
		set logfile [job_logfile compress_$method [pwd] $cmdline {*}[versions $method]]
	}
	foreach file $args {
		if {$outputfile eq ""} {
			set target [gzroot $file].$method
		} else {
			set target $outputfile
		}
		if {[job_distribute] eq "0" || $file eq "-"} {
			compress_$method $file $target $index $keep $threads $compressionlevel $blocksize
		} else {
			job $method-$file -deps {$file} -targets {$target} -vars {method file keep index threads compressionlevel blocksize} -code {
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
