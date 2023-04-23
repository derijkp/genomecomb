proc sam_merge_job {args} {
	upvar job_logdir job_logdir
	set threads 1
	set force 0
	set deletesams 0
	set optional 0
	set index 1
	set sortopts {}
	cg_options sam_merge args {
		-index {
			set index $value
		}
		-deletesams {
			set deletesams $value
		}
		-threads {
			set threads $value
		}
		-force {
			set force $value
		}
		-sort {
			if {$value eq "name"} {set sortopts -n} else {set sortopts {}}
		}
		-optional {
			set optional $value
		}
		-name {
			set name $value
		}
		-m - -maxopenfiles {
			maxopenfiles $value
		}
	} {outputfile samfile} 1 ... {
		merge sam files using samtools merge
	}
	set samfiles [list $samfile {*}$args]
	if {![info exists name]} {
		set name [job_relfile2name sam_merge- $outputfile]
	}
	set outputfile [file_absolute $outputfile]
	if {![jobtargetexists -checkdepexists 1 $outputfile $samfiles]} {
		set workdir [shadow_workdir $outputfile]/merge
		file delete -force $workdir
		file mkdir $workdir
		set_job_logdir [file dir $outputfile]/log_jobs
		set maxfiles [maxopenfiles]
		if {$maxfiles < 2} {set maxfiles 2}
		set len [llength $samfiles]
		if {$len <= $maxfiles} {
			set target $outputfile
			job $name -optional $optional -force $force -deps $samfiles -targets {
				$target
			} -vars {threads sortopts} -code {
				exec samtools merge {*}$sortopts -t $threads $target.temp {*}$deps
				file rename -force -- $target.temp $target
			}
		} else {
			catch {file delete {*}[glob -nocomplain $workdir/paste.temp*]}
			set todo $samfiles
			set delete 0
			set num 1
			while 1 {
				if {$len <= $maxfiles} {
					set target $outputfile
					job paste-$name -optional $optional -force $force -cores $threads \
					-deps $todo -targets {
						$target
					} -vars {delete workdir threads sortopts} -code {
						exec samtools merge {*}$sortopts -t $threads $target.temp {*}$deps
						file rename -force -- $target.temp $target
						if {$delete} {file delete {*}$deps}
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
					job ${name}-$num -optional $optional -force $force \
					-deps $deps -targets {
						$target
					} -vars {delete threads sortopts} -code {
						# puts [list ../bin/tsv_paste {*}$deps]
						if {[llength $deps] > 1} {
							exec samtools merge {*}$sortopts -t $threads $target.temp {*}$deps
							if {$delete} {file delete {*}$deps}
						} elseif {!$delete} {
							mklink $dep $target.temp
						} else {
							file rename $dep $target.temp
						}
						file rename -force -- $target.temp $target
					}
				}
				set delete 1
				set todo $newtodo
				set len [llength $todo]
			}
		}
	}
	if {$deletesams} {	
		job $name-deletesams -optional $optional -force $force \
		-deps [list $outputfile $samfiles] -vars {samfiles} -rmtargets $samfiles -code {
			file delete {*}$samfiles
		}
		
	}
	if {$index} {	
		bam_index_job -optional $optional -force $force -threads $threads $outputfile
	}
}

proc cg_sam_merge {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [sam_merge_job {*}$args]
	job_wait
	return $result
}
