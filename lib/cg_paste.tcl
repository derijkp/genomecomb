proc tsv_paste_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg paste {*}$args]"
	set forcepaste 0
	set endcommand {}
	set optional 0
	cg_options paste args {
		-o - -outputfile {
			set outputfile $value
		}
		-m - -maxopenfiles {
			maxopenfiles $value
		}
		-forcepaste {
			set forcepaste $value
		}
		-endcommand {
			set endcommand $value
		}
		-optional {
			set optional $value
		}
	} {} 0
	if {![info exists outputfile]} {
		# puts [list ../bin/tsv_paste {*}$args]
		exec tsv_paste {*}$args >@ stdout 2>@ stderr
		return
	}
	set files $args
	# putsvars outputfile files forcepaste endcommand
	set outputfile [file_absolute $outputfile]
	set outputdir [file dir $outputfile]
	job_logfile $outputdir/paste_[file tail $outputfile] $outputdir $cmdline
	set maxfiles [maxopenfiles]
	set len [llength $files]
	if {$len <= $maxfiles} {
		set target $outputfile
		job [job_relfile2name paste- $outputfile] -optional $optional -force $forcepaste -deps $files -targets {$target} -vars {endcommand} -code {
			set compress [compresspipe $target]
			# puts [list ../bin/tsv_paste {*}$deps]
			set temp [filetemp_ext $target]
			exec tsv_paste {*}$deps {*}$compress > $temp 2>@ stderr
			file rename -force -- $temp $target
			if {$compress ne ""} {cg_zindex $target}
			if {$endcommand ne ""} {eval $endcommand}
			
		}
		return
	}
	set workdir [workdir $outputfile]/paste
	file delete -force $workdir
	file mkdir $workdir
	catch {file delete {*}[glob -nocomplain $workdir/paste.temp*]}
	set todo $files
	set delete 0
	set num 1
	while 1 {
		if {$len <= $maxfiles} {
			set target $outputfile
			job [job_relfile2name paste- $outputfile] -optional $optional -force $forcepaste \
			-deps $todo -targets {
				$target
			} -vars {
				endcommand delete workdir
			} -code {
				set compress [compresspipe $target]
				# puts [list ../bin/tsv_paste {*}$deps]
				set temp [filetemp_ext $target]
				exec tsv_paste {*}$deps {*}$compress > $temp 2>@ stderr
				file rename -force -- $temp $target
				if {$compress ne ""} {cg_zindex $target}
				if {$delete} {file delete {*}$deps}
				if {$endcommand ne ""} {eval $endcommand}
				file delete -force $workdir
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
			job [job_relfile2name paste- $target] -optional $optional -deps $deps -force $forcepaste -targets {$target} -vars {delete} -code {
				# puts [list ../bin/tsv_paste {*}$deps]
				if {[llength $deps] > 1} {
					exec tsv_paste {*}$deps {*}[compresspipe $target 1] > $target.temp 2>@ stderr
					if {$delete} {file delete {*}$deps}
				} elseif {[file extension $dep] ne ".zst"} {
					exec {*}[compresscmd $target 1 1] $dep > $target.temp
				} elseif {!$delete} {
					mklink $dep $target.temp
				} else {
					file rename -- $dep $target.temp
				}
				file rename -force -- $target.temp $target
			}
		}
		set delete 1
		set todo $newtodo
		set len [llength $todo]
	}
}

proc cg_paste {args} {
	set args [job_init {*}$args]
	tsv_paste_job {*}$args
	job_wait
}

