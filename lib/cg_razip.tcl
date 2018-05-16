proc cg_razip_job args {
	upvar job_logdir job_logdir
	set pos 0
	set keep 0
	set cmdline [list cg razip {*}$args]
	cg_options razip args {
		-k - -keep {
			set keep 1
		}
	}
	set logfile [job_logfile razip [pwd] $cmdline {*}[versions razip]]
	foreach file $args {
		set target [gzroot $file].rz
		job razip-$file -deps {$file} -targets {$target} -vars {keep} -code {
			set file $dep
			set ext [file extension $file]
			switch $ext {
				.gz {
					putslog "razip $file"
					set result [file root $file].rz
					exec gunzip -d -c $file > $result.temp2
					exec razip -c $result.temp2 > $result.temp
					file delete $result.temp2
					file rename -force $result.temp $result
					if {!$keep} {file delete $file}
				}
				.rz {
					putslog "$file already razip"
				}
				.lz4 {
					putslog "razip $file"
					set result [file root $file].rz
					exec lz4c -q -d $file > $result.temp2
					exec razip -c $result.temp2 > $result.temp
					file delete $result.temp2
					file rename -force $result.temp $result
					if {!$keep} {file delete $file}
				}
				.bz2 {
					putslog "razip $file"
					set result [file root $file].rz
					exec bzcat $file > $result.temp2
					exec razip -c $result.temp2 > $result.temp
					file delete $result.temp2
					file rename -force $result.temp $result
					if {!$keep} {file delete $file}
				}
				default {
					putslog "razip $file"
					exec razip -c $file > $file.rz.temp
					file rename -force $file.rz.temp $file.rz
					if {!$keep} {file delete $file}
				}
			}
		}
	}
}

proc cg_razip args {
	set args [job_init {*}$args]
	cg_razip_job {*}$args
	job_wait
}
