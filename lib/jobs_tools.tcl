proc cleanup_job {name rmtargets args} {
	upvar job_logdir job_logdir
	set todo 0
	foreach temp $rmtargets {
		if {[jobfileexists $temp]} {set todo 1}
		set analysisinfo [gzroot $temp].analysisinfo
		if {[jobfileexists $analysisinfo]} {
			set todo 1
			lappend rmtargets $analysisinfo
		}
	}
	if {!$todo} return
	set num 1
	foreach deps $args {
		job cleanup-$name-deps$num -optional 1 -deps [list {*}$deps] -vars {rmtargets} \
		-rmtargets $rmtargets -code {
			foreach file $rmtargets {
				if {[file isdir $file]} {
					if {![llength [glob -nocomplain $file/*]]} {
						catch {file delete -force $file}
					}
				} else {
					catch {file delete -force $file}
				}
			}
		}
		incr num
	}
}

proc job_optdeps {deps} {
	return \([join $deps \)\ \(]\)
}
