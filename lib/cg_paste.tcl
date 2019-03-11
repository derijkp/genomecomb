proc tsv_paste_job {outputfile files args} {
	set forcepaste 0
	set endcommand {}
	set optional 0
	foreach {k v} $args {
		switch $k {
			-forcepaste {set forcepaste $v}
			-endcommand {set endcommand $v}
			-optional {set optional $v}
			default {error "Unkown option $k"}
		} 
	}
	# putsvars outputfile files forcepaste endcommand
	set outputfile [file_absolute $outputfile]
	set workdir [gzroot $outputfile].index/paste
	file delete -force $workdir
	file mkdir $workdir
	job_logdir $workdir/log_jobs
	set maxfiles [maxopenfiles]
	if {$maxfiles < 2} {set maxfiles 2}
	set len [llength $files]
	if {$len <= $maxfiles} {
		set target $outputfile
		job paste-[file tail $outputfile] -optional $optional -force $forcepaste -deps $files -targets {$target} -vars {endcommand} -code {
			set compress [compresspipe $target]
			# puts [list ../bin/tsv_paste {*}$deps]
			set temp [filetemp_ext $target]
			exec tsv_paste {*}$deps {*}$compress > $temp 2>@ stderr
			file rename -force $temp $target
			if {$compress ne ""} {cg_zindex $target}
			if {$endcommand ne ""} {eval $endcommand}
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
			job paste-[file tail $outputfile] -optional $optional -force $forcepaste -deps $todo -targets {$target} -vars {endcommand delete workdir} -code {
				set compress [compresspipe $target]
				# puts [list ../bin/tsv_paste {*}$deps]
				set temp [filetemp_ext $target]
				exec tsv_paste {*}$deps {*}$compress > $temp 2>@ stderr
				file rename -force $temp $target
				if {$compress ne ""} {cg_zindex $target}
				if {$delete} {file delete {*}$deps}
				if {$endcommand ne ""} {eval $endcommand}
			}
			break
		}
		set pos 0
		set newtodo {}
		while {$pos < $len} {
			set target $workdir/paste.temp$num.zst
			incr num
			lappend newtodo $target
			set deps [lrange $todo $pos [expr {$pos+$maxfiles-1}]]
			incr pos $maxfiles
			job paste-[file tail $target] -optional $optional -deps $deps -force $forcepaste -targets {$target} -vars {delete} -code {
				# puts [list ../bin/tsv_paste {*}$deps]
				if {[llength $deps] > 1} {
					exec tsv_paste {*}$deps {*}[compresspipe $target 1] > $target.temp 2>@ stderr
					if {$delete} {file delete {*}$deps}
				} elseif {[file extension $dep] ne ".zst"} {
					exec {*}[compresscmd $target 1 1] $dep > $target.temp
				} elseif {!$delete} {
					mklink $dep $target.temp
				} else {
					file rename $dep $target.temp
				}
				file rename -force $target.temp $target
			}
		}
		set delete 1
		set todo $newtodo
		set len [llength $todo]
	}
}

proc cg_paste {args} {
	set args [job_init {*}$args]
	cg_options paste args {
		-o - -outputfile {
			set outputfile $value
		}
		-m - -maxopenfiles {
			set ::maxopenfiles $value
		}
	} {} 1
	if {[info exists outputfile]} {
		tsv_paste_job $outputfile $args
	} else {
		# puts [list ../bin/tsv_paste {*}$args]
		exec tsv_paste {*}$args >@ stdout 2>@ stderr
	}
}

