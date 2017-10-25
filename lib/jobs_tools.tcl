proc cleanup_job {name rmtargets args} {
	upvar job_logdir job_logdir
	set todo 0
	foreach temp $rmtargets {
		if {[jobfileexists $temp]} {set todo 1; break}
	}
	if {!$todo} return
	set num 1
	foreach deps $args {
		job cleanup-$name-deps$num -optional 1 -deps [list {*}$deps] -vars {rmtargets} \
			-rmtargets $rmtargets -code {
			foreach file $rmtargets {
				catch {file delete $file}
				catch {file delete [gzroot $file].analysisinfo}
			}
		}
		incr num
	}
}

proc job_optdeps {deps} {
	return \([join $deps \)\ \(]\)
}
