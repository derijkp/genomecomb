# profiler start
# # code
# profiler cp1
# # code
# profiler cp2
# # code
# profiler end

proc profiler {args} {
	global db_t
	global db_c
	array set db_t {}
	array set db_c {}
	global last
	set id [lindex $args 0]
	switch -- $id {
		start {
			set last [clock microseconds]
		}
		end {
			set file [lindex $args 1]
			if {$file eq ""} {set file stdout}
			set k [array names db_t]
			puts [format {%-12s %s} {checkpoint:} {avgtime:}]
			foreach ik $k {
				puts $file [format "%s\t%.1f\t%10d\t%.1f" $ik [expr {1.0*$db_t($ik)/$db_c($ik)}] $db_c($ik) $db_t($ik)]
			}
			array unset db_t
			array unset db_c
		}
		default {
			set delta [expr {[clock microseconds]-$last}]
			set last [clock microseconds]
			if {[info exists db_t($id)]} {incr db_t($id) $delta} {set db_t($id) $delta}
			if {[info exists db_c($id)]} {incr db_c($id) 1     } {set db_c($id) 1     }
		}
	}
}
