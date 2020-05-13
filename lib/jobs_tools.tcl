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
		if {[jobfileexists $analysisinfo]} {
			set todo 1
			lappend rmtargets $analysisinfo
		}
		set indexfile [index_file $temp]
		if {[jobfileexists $indexfile]} {
			set todo 1
			lappend rmtargets $indexfile
		}
	}
	if {!$todo} return
	set num 1
	foreach deps $args {
		set deps [list_remove $deps {}]
		job cleanup-$name-deps$num -optional 1 -deps [list {*}$deps] -rmtargets $rmtargets -vars {
			rmtargets forcedirs delassociated
		} -code {
			foreach file $rmtargets {
				if {!$forcedirs && [file isdir $file]} {
					if {![llength [glob -nocomplain $file/*]]} {
						catch {file delete -force $file}
					}
				} else {
					catch {file delete -force $file}
				}
				if {$delassociated} {
					set indexfile [index_file $file]
					if {[file exists $indexfile]} {
						file delete -force $indexfile
					}
					set analysisinfofile [analysisinfo_file $file]
					if {[file exists $analysisinfofile]} {file delete -force $analysisinfofile}
					if {[file exists [gzroot $file].index] && ($forcedirs || $delassociated > 1)} {
						file delete -force [gzroot $file].index
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
