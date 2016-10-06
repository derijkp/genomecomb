proc tsv_paste_job {outputfile files {forcepaste 0} {endcmd {}}} {
	set outputfile [file_absolute $outputfile]
	set workdir $outputfile.index/paste
	file delete -force $workdir
	file mkdir $workdir
	job_logdir $workdir/log_jobs
	set maxfiles [maxopenfiles]
	if {$maxfiles < 2} {set maxfiles 2}
	set len [llength $files]
	if {$len <= $maxfiles} {
		set target $outputfile
		job paste-[file tail $outputfile] -force $forcepaste -deps $files -targets {$target} -vars {endcmd} -code {
			# puts [list ../bin/tsv_paste {*}$deps]
			exec tsv_paste {*}$deps > $target.temp 2>@ stderr
			file rename -force $target.temp $target
			if {$endcmd ne ""} {eval $endcmd}
		}
		return
	}
	catch {file delete {*}[glob -nocomplain $workdir/paste.temp*]}
	set todo $files
	set delete 0
	set num 1
	while 1 {
		if {$len <= $maxfiles} {
			set target $outputfile
			job paste-[file tail $outputfile] -force $forcepaste -deps $todo -targets {$target} -vars {endcmd delete workdir} -code {
				# puts [list ../bin/tsv_paste {*}$deps]
				exec tsv_paste {*}$deps > $target.temp 2>@ stderr
				file rename -force $target.temp $target
				if {$delete} {file delete {*}$deps $workdir}
				if {$endcmd ne ""} {eval $endcmd}
			}
			break
		}
		set pos 0
		set newtodo {}
		while {$pos < $len} {
			set target $workdir/paste.temp$num
			incr num
			lappend newtodo $target
			set deps [lrange $todo $pos [expr {$pos+$maxfiles-1}]]
			incr pos $maxfiles
			job paste-[file tail $target] -deps $deps -targets {$target} -vars {delete} -code {
				# puts [list ../bin/tsv_paste {*}$deps]
				if {[llength $deps] > 1} {
					exec tsv_paste {*}$deps > $target.temp 2>@ stderr
				} elseif {!$delete} {
					mklink $dep $target.temp
				} else {
					file rename $dep $target.temp
				}
				file rename -force $target.temp $target
				if {$delete} {file delete {*}$deps}
			}
			
		}
		set delete 1
		set todo $newtodo
		set len [llength $todo]
	}
}

proc cg_paste {args} {
	set args [job_init -silent {*}$args]
	cg_options paste args {
		-o - --outputfile {
			set outputfile $value
		}
		-m - --maxopenfiles {
			set ::maxopenfiles $value
		}
	} 1
	if {[info exists outputfile]} {
		tsv_paste_job $outputfile $args
	} else {
		# puts [list ../bin/tsv_paste {*}$args]
		exec tsv_paste {*}$args >@ stdout 2>@ stderr
	}
}

