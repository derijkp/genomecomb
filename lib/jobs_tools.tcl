proc cleanup_job {args} {
	upvar job_logdir job_logdir
	set forcedirs 0
	set delassociated 0
	cg_options cleanup_job args {
		-forcedirs {
			set forcedirs $value
		}
		-delassociated {
			set delassociated $value
		}
	} {name rmtargets} 2
	set rmtargets [list_remove $rmtargets {}]
	set todo 0
	foreach temp $rmtargets {
		if {[jobfileexists $temp]} {set todo 1}
		set analysisinfo [analysisinfo_file $temp]
		if {[jobfileexists -checkcompressed 1 $analysisinfo]} {
			set todo 1
			lappend rmtargets $analysisinfo
		}
		set indexfile [index_file $temp]
		if {[jobfileexists -checkcompressed 1 $indexfile]} {
			set todo 1
			lappend rmtargets $indexfile
		}
	}
	if {!$todo} return
	set num 1
	foreach deps $args {
		set deps [list_remove $deps {}]
		set deps \([join $deps "\)\ \("]\)
		job cleanup-$name-deps$num -optional 1 -checkcompressed 1 -deps [list {*}$deps] -rmtargets $rmtargets -vars {
			rmtargets forcedirs delassociated
		} -code {
			foreach file $rmtargets {
				if {$forcedirs} {
					shadow_delete $file
				} else {
					job_delete_ifempty $file
				}
				if {$delassociated} {
					set indexfile [index_file $file]
					if {[file exists $indexfile]} {
						shadow_delete $indexfile
					}
					set analysisinfofile [analysisinfo_file $file]
					if {[file exists $analysisinfofile]} {shadow_delete $analysisinfofile}
					if {[file exists [gzroot $file].index] && ($forcedirs || $delassociated > 1)} {
						shadow_delete [gzroot $file].index
					}
				}
			}
		}
		incr num
	}
}

proc job_optdeps {deps} {
	return \([join $deps \)\ \(]\)
}

proc job_reset_olds {dir {exclude log_jobs}} {
	# putsvars dir
	foreach file [glob -nocomplain $dir/*] {
		if {[file tail $file] in $exclude} {
			puts "excluded $file"
		} elseif {[regexp \\.old$ $file]} {
			if {![file exists [file root $file]]} {
				puts [list mv $file [file root $file]]
				file rename $file [file root $file]
			} else {
				puts "skipped $file (exists without .old)"
			}
		} elseif {[file isdir $file]} {
			job_reset_olds $file $exclude
		}
	}
}

