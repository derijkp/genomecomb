proc cleanup_job {name rmtargets args} {
	upvar job_logdir job_logdir
	set num 1
	foreach deps $args {
		job cleanup-$name-deps$num -deps [list {*}$deps {*}[job_optdeps $rmtargets]] -vars {rmtargets} \
			-rmtargets $rmtargets -code {
			foreach file $rmtargets {
				catch {file delete $file}
			}
		}
		incr num
	}
}

proc job_optdeps {deps} {
	return \([join $deps \)\ \(]\)
}
