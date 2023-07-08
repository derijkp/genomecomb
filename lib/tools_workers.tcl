proc worker_start {{num 1}} {
	global workers workersa
	set result {}
	for {set i 0} {$i < $num} {incr i} {
		set w [open "|cgworker" r+]
		fconfigure $w -buffering line
		lappend workers $w
		lappend result $w
		set workersa($w) e
	}
	return $result
}

proc worker_close {{ws {}}} {
	global workers
	if {$ws eq {}} {set ws $workers}
	foreach w $ws {
		catch {close $w}
		set workers [list_remove $workers $w]
	}
	return $workers
}

proc worker_exec {w cmd} {
	global workersa
	if {$workersa($w) eq "r"} {
		error "worker $w has a result, \"call worker_get $w\" before calling new worker_exec"
	}
	set send [split $cmd \n]
	puts $w $send
	# puts $w $cmd
	set workersa($w) r
}

proc worker_get {w} {
	global workersa
	if {$workersa($w) eq "e"} {
		error "no result available yet for worker $w: call worker_exec before calling worker_get"
	}
	set result [gets $w]
	set workersa($w) e
	set error [lindex $result 0]
	set result [join [lrange $result 1 end] \n]
	# set result [lindex $result 1]
	if {$error} {
		error "worker $w error: $result"
	}
	return $result
}

proc workers {{num {}}} {
	global workers
	if {![info exists workers]} {set workers {}}
	if {$num ne ""} {
		set len [llength $workers]
		if {$num > $len} {
			worker_start [expr {$num-$len}]
		}
		return [lrange $workers 0 [expr {$num-1}]]
	} else {
		return $workers
	}
}

