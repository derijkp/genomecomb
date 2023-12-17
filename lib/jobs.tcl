proc job.file {type jobfile} {
	mkdir [file dir $jobfile]/$type
	return [file dir $jobfile]/$type/[file tail $jobfile].$type
}

proc job_file_or_link_exists {file} {
	if {[file exists $file]} {return 1}
	if {[catch {file link $file}]} {return 0} else {return 1}
}

proc job_file_mtime {file} {
	file lstat $file a
	return $a(mtime)
}

proc job_distribute {{type {}}} {
	global cgjob
	if {$type eq ""} {
		return $cgjob(distribute)
	}
	set cgjob(distribute) $type
	if {[isint $cgjob(distribute)]} {
		if {$cgjob(distribute) <= 1} {
			set target direct
		} else {
			set target distr
		}
	} else {
		if {$type eq "distr"} {
			set options [list_remove [list_regsub ^job_process_init_ [command_list job_process_init_*] {}] distr]
			error "$type is not a valid option for job distribution, must be an integer (numbner of cores/threads) or one of: sge slurm direct status"
		}
		set target $type
	}
	if {[info commands job_process_init_${target}] eq ""} {auto_load job_process_init_${target}}
#	interp alias {} job_process {} job_process_$target
#	interp alias {} job_wait {} job_wait_${target}
	if {![command_exists job_process_init_${target}]} {
		set options [list_remove [list_regsub ^job_process_init_ [command_list job_process_init_*] {}] distr]
		error "$target is not a valid option for job distribution, must be an integer (numbner of cores/threads) or one of: sge slurm direct status"
	}
	job_process_init_${target}
}

proc job_force_get {{def {}}} {
	global cgjob
	if {$def eq ""} {return $cgjob(force)}
	if {![info exists cgjob(forceset)]} {return $def} else {return $::cgjob(force)}
}

proc job_args {jobargs} {
	global cgjob
	upvar job_logdir job_logdir
	if {![info exists cgjob(distribute)]} {
		set cgjob(distribute) 0
	}
	if {![info exists cgjob(force)]} {
		set cgjob(force) 0
	}
	if {![info exists cgjob(cleanup)]} {
		set cgjob(cleanup) success
	}
	if {![info exists cgjob(removeold)]} {
		set cgjob(removeold) 0
	}
	if {![info exists cgjob(silent)]} {
		if {[get ::verbose 0] >= 1} {
			set cgjob(silent) 0
		} else {
			set cgjob(silent) 1
		}
	}
	if {![info exists cgjob(dry)]} {
		set cgjob(dry) 0
	}
	if {![info exists cgjob(debug)]} {
		set cgjob(debug) 0
	}
	if {![info exists cgjob(nosubmit)]} {
		set cgjob(nosubmit) 0
	}
	if {![info exists cgjob(resubmit)]} {
		set cgjob(resubmit) 1
	}
	if {![info exists cgjob(skipjoberrors)]} {
		set cgjob(skipjoberrors) 0
	}
	if {![info exists cgjob(priority)]} {
		set cgjob(priority) 0
	}
	if {![info exists cgjob(dmem)]} {
		set cgjob(dmem) {}
	}
	if {![info exists cgjob(dtime)]} {
		set cgjob(dtime) {}
	}
#	if {![info exists cgjob(dqueue)]} {
#		set cgjob(dqueue) all.q
#	}
	if {![llength $jobargs]} {return {}}
	set newargs {}
	set pos 0
	while {$pos < [llength $jobargs]} {
		set key [lindex $jobargs $pos]
		incr pos
		switch -- $key {
			-f - -force - --force {
				set cgjob(force) [lindex $jobargs $pos]
				set cgjob(forceset) 1
				incr pos
			}
			-d - -distribute - --distribute {
				job_distribute [lindex $jobargs $pos]
				incr pos
				set cgjob(hasargs) 1
			}
			-dpriority {
				set cgjob(priority) [lindex $jobargs $pos]
				incr pos
			}
			-dqueue {
				set cgjob(dqueue) [lindex $jobargs $pos]
				incr pos
			}
			-dmem {
				set cgjob(dmem) [lindex $jobargs $pos]
				incr pos
			}
			-dtime {
				set cgjob(dtime) [lindex $jobargs $pos]
				incr pos
			}
			-dcleanup {
				set value [lindex $jobargs $pos]
				incr pos
				if {$value ni {success never allways}} {error "$value not a valid option for -dcleanup, should be one of: success, never, allways"}
				set cgjob(cleanup) $value
			}
			-dremoveold {
				set value [lindex $jobargs $pos]
				incr pos
				if {$value ni {0 1}} {error "$value not a valid option for -dremoveold, should be one of: 0 1"}
				set cgjob(removeold) $value
			}
			-silent - --silent {
				set val [lindex $jobargs $pos]
				if {[inlist {0 1} $val]} {
					if {$val == 0} {logverbose 0}
					set cgjob(silent) $val
					incr pos
				} else {
					set cgjob(silent) 1
				}
			}
			-runcmd - --runcmd - -runcommand {
				set cgjob(runcmd) [lindex $jobargs $pos]
				incr pos
			}
			-dry {
				set val [lindex $jobargs $pos]
				set cgjob(dry) $val
				incr pos
			}
			-debug - --debug {
				set val [lindex $jobargs $pos]
				set cgjob(debug) $val
			}
			-dnosubmit {
				set cgjob(nosubmit) [lindex $jobargs $pos]
				incr pos
			}
			-noresubmit - --noresubmit {
				set cgjob(resubmit) 0
			}
			-skipjoberrors - --skipjoberrors {
				set val [lindex $jobargs $pos]
				if {[inlist {0 1} $val]} {
					set cgjob(skipjoberrors) $val
					incr pos
				} else {
					set cgjob(skipjoberrors) 1
				}
			}
			-logfile {
				set logfile [lindex $jobargs $pos]
				incr pos
				job_logfile $logfile
				set_job_logdir [file_absolute [file dir $logfile]/log_jobs]
			}
			-v - -verbose - --verbose {
				set value [lindex $jobargs $pos]
				incr pos
				if {![isint $value]} {error "$value is not a number, only numbers are accepted as value for -v (--verbose)"}
				logverbose $value
			}
			-stack - --stack {
				set value [lindex $jobargs $pos]
				incr pos
				if {$value in {0 1}} {
					set ::stacktraceonerror $value
				} else {
					error "$value is not 0 or 1, the only accepted values for --stack"
				}
			}
			-- {
				lappend newargs --
				break
			}
			default {
				lappend newargs $key
			}
		}
	}
	lappend newargs {*}[lrange $jobargs $pos end]
}

# job_expandvars fills in all variables (from one level up)
# variables prefixed with _ are considered a list, and will be "expanded"
# e.g.
# set a A ; set b {1 2}
# job_expandvars {test-$a$_b}
# test-A1 test-A2
proc job_expandvars {string {level 0}} {
	# find which vars are by putting a subst of the string in a proc and checking the error
	# add set variable value as needed
	# expanding _ prefixed vars may cause multiple results, keep in resultlist, combining 
	# different _ vars as needed
	incr level
	set code "proc job_testvars \{string\} \{"
	set list [list $code 0 {}]
	set resultlist {}
	while 1 {
		set alldone 1
		foreach {code finished result} $list {
			if {$finished} {
				lappend resultlist $code $finished $result
				continue
			}
			set alldone 0
			set cmd "$code\nsubst -nobackslashes -nocommands \$string\n\}"
			eval $cmd
			if {![catch {job_testvars $string} result]} {
				lappend resultlist $code 1 $result
				continue
			}
			if {![regexp {can't read "(.*)": no such variable} $result temp var]} {
				putslog "cannot expand $string: ERROR: $result"
				error $result $::errorInfo
			}
			if {[string index $var 0] eq "_"} {
				set var [string range $var 1 end]
				set expand 1
			} else {
				set expand 0
			}
			if {![uplevel $level info exists $var]} {
				error "Error: can't read \"$var\": no such variable $var"
			}
			set value [uplevel $level get $var]
			if {$expand} {
				foreach value $value {
					lappend resultlist "$code\n[list set _$var $value]" 0 {}
				}
			} else {
				append code "\n[list set $var $value]"
				lappend resultlist $code 0 {}
			}
		}
		if {$alldone} break
		set list $resultlist
		set resultlist {}
	}
	set result {}
	foreach {code finished res} $resultlist { 
		lappend result $res
	}
	return $result
}

proc job_expandvarslist {list {level 0}} {
	set result {}
	incr level
	set list [string_change $list [list \\ \\\\]]
	foreach string $list {
		lappend result {*}[job_expandvars $string $level]
	}
	return [list_remove $result {}]
}

proc jobglob {args} {
	set checkcompressed 1
	set pos 0
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-checkcompressed {
				incr pos
				set checkcompressed [lindex $args $pos]
			}
			-- {
				incr pos
				break
			}
			default break
		}
		incr pos
	}
	set args [lrange $args $pos end]
	set resultfiles {}
	set ids {}
	set time now
	foreach pattern $args {
		set files [job_finddep $pattern ids time timefile $checkcompressed regexppattern]
		if {$regexppattern} {
			set pattern [string range $pattern 1 end-1]
		}
		if {[file pathtype $pattern] eq "relative"} {set relative 1} else {set relative 0}
		foreach file $files {
			if {$relative} {regsub ^[pwd]/ $file {} file}
			lappend resultfiles $file
		}
	}
	list_remdup $resultfiles
}

proc jobgzfiles {args} {
	foreach filename $args {
		set list [jobglob $filename $filename.zst $filename.lz4 $filename.rz $filename.bgz $filename.gz $filename.bz2]
		foreach file $list {
			set root [gzroot $file]
			if {[info exists a($root)]} continue
			set a($root) $file
		}
	}
	set result {}
	foreach file [array names a] {
		lappend result $a($file)
	}
	return $result
}

proc jobgzfile {args} {
	set list [jobgzfiles {*}$args]
	if {[llength $list]} {
		return [lindex $list 0]
	} else {
		return [lindex $args 0]
	}
}

proc jobglob1 {args} {
	lindex [jobglob {*}$args] 0
}

proc jobfileexists {args} {
	set checkcompressed 1
	set pos 0
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-checkcompressed {
				incr pos
				set checkcompressed [lindex $args $pos]
			}
			-- {
				incr pos
				break
			}
			default break
		}
		incr pos
	}
	set args [lrange $args $pos end]
	set files {}
	set ids {}
	set time now
	foreach pattern $args {
		if {[file exists $pattern]} {
		} elseif {![llength [job_finddep $pattern ids time timefile $checkcompressed]]} {
			return 0
		}
	}
	return 1
}

# jobtargetexists ?options? target deps
# returns 1 if target exists, and non of the dependencies (args) is newer than target (or being made)
# if -checkdepexists is 0 (default), dependencies are only checked if they exist
proc jobtargetexists {args} {
	set checkcompressed 1
	set checkdepexists 0
	set pos 0
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-checkcompressed {
				incr pos
				set checkcompressed [lindex $args $pos]
			}
			-checkdepexists {
				incr pos
				set checkdepexists [lindex $args $pos]
			}
			-- {
				incr pos
				break
			}
			default break
		}
		incr pos
	}
	set args [lrange $args $pos end]
	foreach {targets deps} $args break
	set targettime new
	foreach target $targets {
		if {![file exists $target]} {return 0}
		if {$targettime eq "new" || [job_file_mtime $target] < $targettime} {
			set targettime [job_file_mtime $target]
		}
	}
	foreach pattern $deps {
		set time 0
		set files [job_finddep $pattern ids time timefile $checkcompressed]
		if {($checkdepexists && ![llength $files]) || $time in "now force" || $time > $targettime} {
#			job_log [job_relfile2name targetexists- $target] "one of targets older than dep $timefile (renaming to .old): $targets"
#			foreach target $targets {
#				job_to_old $target
#			}
			return 0
		}
	}
	return 1
}

proc glob2regexp {pattern} {
	set regexp {}
	set escape 0
	set pos 0
	set len [string length $pattern]
	while {$pos < $len} {
		set c [string index $pattern $pos]
		incr pos
		if {$escape} {
			append regexp $c
			set escape 0
		} elseif {$c eq "\\"} {
			append regexp $c
			set escape 1
		} elseif {$c eq "*"} {
			append regexp {(.*)}
		} elseif {$c eq "?"} {
			append regexp {(.)}
		} elseif {$c eq "\["} {
			append regexp "\(\["
		} elseif {$c eq "\]"} {
			append regexp "\]\)"
		} elseif {$c eq "\{"} {
			append regexp "\("
			while {$pos < $len} {
				set c [string index $pattern $pos]
				incr pos
				if {$c eq "\}"} break
				if {$c eq "\,"} {
					append regexp "|"
				} elseif {[regexp {[A-Za-z0-9 ]} $c]} {
					append regexp "$c"
				} else {
					append regexp "\\$c"
				}
			}
			append regexp "\)"
		} elseif {[regexp {[A-Za-z0-9 ]} $c]} {
			append regexp "$c"
		} else {
			append regexp "\\$c"
		}
	}
	return $regexp
}

if 0 {
	set p [glob2regexp ab*.cd{ab,cd}?f]
	regexp $p abxx.cdcdef
	regexp $p abxx.cdabef
	regexp $p abxxcdcdef
}

proc cgjob_files {cgjob_idVar pattern {checkcompressed 0}} {
	upvar $cgjob_idVar cgjob_id
	set filelist {}
	regsub -all {([^A-Za-z0-9*])} $pattern {\\\1} rpattern
	if {$checkcompressed} {
		set list [array names cgjob_id $pattern*]
		set regexppattern ^[string_change $rpattern {* [^/]*}](\.zst|\.gz|\.rz|\.lz4|\.bgz|\.bz2|)\$
	} else {
		set list [array names cgjob_id $pattern]
		set regexppattern ^[string_change $rpattern {* [^/]*}]\$
	}
	foreach file $list {
		if {[regexp $regexppattern $file]} {
			lappend filelist $file
		}
	}
	return $filelist
}

proc job_finddep {pattern idsVar timeVar timefileVar checkcompressed {regexppatternVar {}}} {
	global cgjob_id cgjob_rm
#puts *****************
#putsvars pattern checkcompressed
#puts cgjob_id:[array names cgjob_id]
#puts cgjob_rm:[array names cgjob_rm]
	upvar $idsVar ids
	upvar $timeVar time
	upvar $timefileVar timefile
	if {$regexppatternVar ne ""} {upvar $regexppatternVar regexppattern}
	set regexppattern 0
	if {$pattern eq ""} {return {}}
	if {[string index $pattern 0] eq "^" && [string index $pattern end] eq "\$"} {
		set regexppattern 1
		set pattern [string range $pattern 1 end-1]
		return [job_findregexpdep $pattern ids time timefile $checkcompressed]
	}
	set pattern [file_absolute $pattern]
	unset -nocomplain filesa
	if {$checkcompressed} {
		set list [bsort [gzfiles $pattern]]
	} else {
		set list [bsort [checkfiles $pattern]]
	}
	foreach file $list {
		if {[info exists cgjob_rm($file)]} continue
		if {[gziscompressed $file] && [info exists filesa([file root $file])]} {
			continue
		}
		set filesa($file) 1
		lappend ids {}
		maxfiletime $file time timefile
	}
	set filelist [cgjob_files cgjob_id $pattern $checkcompressed]
	foreach file $filelist {
		if {[info exists cgjob_rm($file)]} continue
		if {[gziscompressed $file] && [info exists filesa([file root $file])]} {
			continue
		}
# do not remove job, we want to keep the dependency chain intact,
# even if one job stops prematurely
#		if {![job_running $cgjob_id($file)]} {
#			unset cgjob_id($file)
#			continue
#		}
		set filesa($file) 1
		lappend ids $cgjob_id($file)
		if {[job_running $cgjob_id($file)]} {
			set time now
			set timefile $file
		}
	}
	return [array names filesa]
}

proc maxfiletime {file timeVar timefileVar} {
	upvar $timeVar time
	upvar $timefileVar timefile
	if {$time eq "now"} {return $time}
	if {$time eq "force"} {return $time}
	set ftime [job_file_mtime $file]
	if {$ftime > $time} {
		set time $ftime
		set timefile $file
	}
}

proc job_findregexpdep {pattern idsVar timeVar timefileVar checkcompressed} {
	global cgjob_id cgjob_rm
	upvar $idsVar ids
	upvar $timeVar time
	upvar $timefileVar timefile
	set pattern [file_absolute $pattern]
	set glob [regexp2glob $pattern]
	unset -nocomplain filesa
	# check file system
	if {$checkcompressed} {
		set list [bsort [gzfiles $glob]]
	} else {
		set list [bsort [checkfiles $glob]]
	}
	foreach file $list {
		if {[info exists cgjob_rm($file)]} continue
		if {[gziscompressed $file] && [info exists filesa([file root $file])]} {
			continue
		}
		if {[regexp ^$pattern\$ $file]} {
			maxfiletime $file time timefile
			set filesa($file) 1
			lappend ids {}
		}
	}
	# check files from running/submitted jobs (cgjob_id)
	if {$checkcompressed} {
		set filelist [gzarraynames cgjob_id [file_absolute $glob]]
	} else {
		set filelist [array names cgjob_id [file_absolute $glob]]
	}
	foreach file $filelist {
		if {[info exists cgjob_rm($file)]} continue
		if {[gziscompressed $file] && [info exists filesa([file root $file])]} {
			continue
		}
		if {![regexp ^[file_absolute $pattern](\.gz|\.rz|\.zst\.lz4|\.bz2|\.bgz)?\$ $file]} continue
# do not remove job, we want to keep the dependency chain intact,
# even if one job stops prematurely
#		if {![job_running $cgjob_id($file)]} {
#			unset cgjob_id($file)
#			continue
#		}
		set filesa($file) 1
		lappend ids $cgjob_id($file)
		if {[job_running $cgjob_id($file)]} {
			set time now
			set timefile $file
		}
	}
	return [array names filesa]
}

# dependencies between braces () are optional (braces must be at start and end of dependency)
# $targetvarsVar will contain the values extracted from () matches in the dep pattern, that will be filled in
# for \1, ... in the target
proc job_finddeps {job deps targetvarsVar targetvarslist idsVar timeVar timefileVar checkcompressed {ftargetvars {}}} {
	upvar $idsVar ids
	if {$targetvarsVar ne ""} {
		upvar $targetvarsVar targetvars
	}
	upvar $timeVar time
	upvar $timefileVar timefile
	set time 0
	set timefile {}
	set ids {}
	set finaldeps {}
	set targetvars {}
	foreach pattern $deps {
		if {[string index $pattern 0] eq "\(" && [string index $pattern end] eq "\)"} {
			set opt 1
			set pattern [string range $pattern 1 end-1]
		} else {
			set opt 0
		}
		set pattern [job_targetreplace $pattern $ftargetvars]
		set files [job_finddep $pattern ids time timefile $checkcompressed regexppattern]
		if ($regexppattern) {
			set pattern [string range $pattern 1 end-1]
		}
		if {![llength $files]} {
			if {!$opt} {
				error "missing dependency \"$pattern\""
			} else {
				job_lognf $job "missing optional dependency \"$pattern\""
				lappend finaldeps {}
			}
			continue
		} else {
			job_lognf $job "dependency ok ($pattern): $files"
		}
		foreach file [list_remdup $files] {
			lappend finaldeps $file
			if {$regexppattern} {
				set targets [lrange [regexp -all -inline ^[file_absolute $pattern]\$ $file] 1 end]
			} else {
				# todo
				set targets [lrange [regexp -all -inline ^[file_absolute [glob2regexp $pattern]]\$ $file] 1 end]
			}
			if {!$targetvarslist} {
				lappend targetvars {*}$targets
			} else {
				# for foreach deps, we must get separate lists for each file
				# in order to be able to insert the correct ftargetvars into the new jobs in the queue
				lappend targetvars $targets
			}
		}
	}
	return $finaldeps
}

# set string {abc\1def\\1ghi\\\1jkl\\\\1s\_a} ; set targetvars {x xx}
proc job_targetreplace {string targetvars {expandVar {}}} {
	set targetvarslen [llength $targetvars]
	set expandlist {}
	if {$expandVar != ""} {
		upvar $expandVar expand
		set expand {}
		set pattern1 {\\+(?:[0-9]+|_)}
		set pattern2 {(\\+)([0-9]+|_)}
	} else {
		set pattern1 {\\+[0-9]+}
		set pattern2 {(\\+)([0-9]+)}
	}
	list_foreach {b e} [list_reverse [regexp -inline -all -indices $pattern1 $string]] {
		regexp $pattern2 [string range $string $b $e] temp pre num
		set len [string length $pre]
		if {![expr {$len%2}]} continue
		if {$num eq "_"} {
			lappend expandlist [string range $string [expr {$e+1}] end]
			set string [string range $string 0 [expr {$b+$len-2}]]
		} else {
			incr num -1
			if {!$targetvarslen} {set value *} else {set value [lindex $targetvars $num]}
			set string [string_replace $string [expr {$b+$len-1}] $e $value]
		}
	}
	if {[llength $expandlist]} {
		lappend expandlist $string
		set expandlist [list_reverse $expandlist]
		set expand {}
		foreach t $targetvars {
			lappend expand [file_absolute [join $expandlist $t]]
		}
		return {}
	} else {
		return $string
	}
}

proc job_targetsreplace {list targetvars} {
	set result {}
	foreach string $list {
		set temp [job_targetreplace $string $targetvars expand]
		if {[llength $expand]} {
			lappend result {*}$expand
		} else {
			lappend result [file_absolute $temp]
		}
	}
	return $result
}

# returns 1 if target is ok (no rerun is needed)
proc job_checktarget {job target skiptarget time timefile checkcompressed {newidsVar {}}} {
	if {$target eq ""} {return 1}
	global cgjob_id
	if {$skiptarget} {set skiptext skip} else {set skiptext ""}
	if {$newidsVar ne ""} {
		upvar $newidsVar newids
	}
	if {$checkcompressed} {
		set files [gzfiles $target]
	} else {
		set files [checkfiles $target]
	}
	if {[llength $files]} {
		if {$time ne "now" && $time ne "force"} {
			foreach file $files {
				if {[job_file_mtime $file] < $time} {
					set time now
					break
				}
			}
		}
		if {$time eq "now"} {
			if {!$skiptarget} {
				job_lognf $job "target older than dep $timefile (renaming to .old): $target"
				foreach file $files {
					job_to_old $file
				}
			} else {
				job_lognf $job "skiptarget older than dep $timefile: $target"
			}
			return 0
		} elseif {$time eq "force"} {
			job_lognf $job "${skiptext}target overwrite (force): $target"
			unset -nocomplain cgjob_id($target)
			return 1
		} else {		
			job_lognf $job "${skiptext}target ok: $target"
			unset -nocomplain cgjob_id($target)
			return 1
		}
	} elseif {[info exists cgjob_id($target)] && $cgjob_id($target) != "q"} {
		job_lognf $job "${skiptext}target already submitted/running (id $cgjob_id($target)): $target "
		set newids $cgjob_id($target)
		return 2
	}
	job_lognf $job "${skiptext}target missing: $target"
	return 0
}

# returns 1 if targets are ok (no rerun is needed)
proc job_checktargets {job targets skiptarget time timefile checkcompressed {runningVar {}}} {
	if {$runningVar ne ""} {
		upvar $runningVar running
		set running {}
	}
	set ok 1
	foreach target $targets {
		set check [job_checktarget $job $target $skiptarget $time $timefile $checkcompressed]
		if {!$check} {
			set ok 0
		}
		if {$check == 2} {
			lappend running $target
		}
	}
	return $ok
}

proc job_logdir {{logdir {}}} {
	upvar job_logdir job_logdir
	if {$logdir eq ""} {
		set logdir [file join [pwd] log_jobs]
	}
	set job_logdir [file_absolute $logdir]
	shadow_mkdir $job_logdir
	set ::cgjob(default_job_logdir) 0
}

proc set_job_logdir {{logdir {}}} {
	upvar job_logdir job_logdir
	if {[info exists job_logdir] && [get ::cgjob(default_job_logdir) 1] == 0} return
	if {$logdir eq ""} {
		set logdir [file join [pwd] log_jobs]
	}
	set job_logdir [file_absolute $logdir]
	shadow_mkdir $job_logdir
	set ::cgjob(default_job_logdir) 0
}

proc job_logname {job_logdir name} {
	set name [string_change $name {/ __ : _ \" _ \' _}]
	set name [shorten_filename $name]
	return [file_absolute $job_logdir/$name]
}

proc job_logclear {job} {
	set ::cgjob(buffer,$job) {}
}

proc job_timestamp {} {
	set now [clock milliseconds]
	set seconds [expr {$now/1000}]
	set milliseconds [expr {$now%1000}]
	return [clock format $seconds -format "%Y-%m-%d %H:%M:%S"].[format %03d $milliseconds]
}

proc job_log {job args} {
	global cgjob
	file mkdir [file dir $job]
	foreach message $args {
		lappend cgjob(buffer,$job) "[job_timestamp]\t$message"
	}
	if {![llength $cgjob(buffer,$job)]} return
	if {![info exists cgjob(f,$job)]} {
		set cgjob(f,$job) [open [job.file log $job] a]
	}
	set f $cgjob(f,$job)
	set log [join $cgjob(buffer,$job) \n]
	puts $f $log
	if {!$cgjob(silent)} {puts stderr $log}
	set cgjob(buffer,$job) {}
	flush $f
}

proc job_getfromlog {job key} {
	set c [file_read [job.file log $job]]
	set line [lindex [regexp -all -inline "$key: \[^\\n\]+" $c] end]
	return [lindex [string range $line [expr {[string length $key]+1}] end] 0]
}

# log to buffer only, only go to output at job_log
# but can be cleared before output using job_logclear 
proc job_lognf {job args} {
	global cgjob
	foreach message $args {
		lappend cgjob(buffer,$job) "[job_timestamp]\t$message"
	}
}

proc job_logclose {job args} {
	global cgjob
	if {![info exists cgjob(f,$job)]} {
		set cgjob(f,$job) [open [job.file log $job] a]
	}
	job_log $job
	set f $cgjob(f,$job)
	catch {close $f}
	unset cgjob(f,$job)
	unset cgjob(buffer,$job)
}

proc job_logfile_set {logfile {dir {}} {cmdline {}} args} {
	global cgjob
	upvar job_logdir job_logdir
	# allways set logdir next to logfile
	set_job_logdir [file dir $logfile]/log_jobs
	if {$dir eq ""} {set dir [file dir $logfile]}
	set time [string_change [timestamp] {" " _ : - . -}]
	set cgjob(logfile) [file_absolute $logfile].$time
	set cgjob(prefix) [file tail $cgjob(logfile)]
	file mkdir [file dir $cgjob(logfile)]
	set cgjob(f_logfile) [open $cgjob(logfile).submitting w]
	puts $cgjob(f_logfile) "\# genomecomb log file"
	set cgjob(basedir) $dir
	if {$dir ne ""} {
		puts $cgjob(f_logfile) "\# basedir: $dir"
	}
	if {$cmdline ne ""} {
		puts $cgjob(f_logfile) "\# cmdline: $cmdline"
	}
	puts $cgjob(f_logfile) "\# version_genomecomb: [version genomecomb]"
	puts $cgjob(f_logfile) "\# distribute: $cgjob(distribute)"
	foreach {key value} $::job_method_info {
		puts $cgjob(f_logfile) "\# $key: $value"
	}
	foreach {key value} $args {
		puts $cgjob(f_logfile) "\# version_$key: $value"
	}
	puts $cgjob(f_logfile) [join {job jobid status submittime starttime endtime duration time_seconds targets msg run cores} \t]
	flush $cgjob(f_logfile)
	set cgjob(totalduration) {0 0}
	set cgjob(status) ok
	set cgjob(starttime) [timestamp]
	return $cgjob(logfile)
}

proc job_logfile {{logfile {}} {dir {}} {cmdline {}} args} {
	global cgjob
	upvar job_logdir job_logdir
	# This will only set the log file the first time it is called, allowing subcommands to set it if called separately
	# but not when called from a larger workflow
	if {$cgjob(logfile) ne ""} {return $cgjob(logfile)}
	if {$dir ne ""} {
		putslog "basedir: $dir"
	}
	if {$cmdline ne ""} {
		putslog "cmdline: $cmdline"
	}
	job_logfile_set $logfile $dir $cmdline {*}$args
	putslog "version_genomecomb: [version genomecomb]"
	putslog "distribute: $cgjob(distribute)"
	putslog "logfile: $cgjob(logfile).*"
	if {$args ne ""} {
		putslog "versions: $args"
	}
}

proc timediff2duration {diff} {
	if {$diff eq ""} {return ""}
	foreach {days miliseconds} $diff break
	set seconds [expr {$days*86400 + $miliseconds/1000.0}]
	set hours [expr {int($seconds/3600)}]
	set minutes [expr {int(($seconds-3600*$hours)/60)}]
	set seconds [expr {$seconds-3600*$hours-60*$minutes}]
	foreach {seconds miliseconds} [split [format %.3f $seconds] .] break
	return $hours:$minutes:$seconds.$miliseconds
}

proc timebetween_induration {starttime endtime} {
	if {[catch {time_scan $endtime} endcode]} {return {}}
	if {[catch {time_scan $starttime} startcode]} {return {}}
	set diff [lmath_calc $endcode - $startcode]
	timediff2duration $diff
}

proc timebetween_inseconds {starttime endtime} {
	if {[catch {time_scan $endtime} endcode]} {return {}}
	if {[catch {time_scan $starttime} startcode]} {return {}}
	set diff [lmath_calc $endcode - $startcode]
	foreach {days miliseconds} $diff break
	format %.3f [expr {$days*86400 + $miliseconds/1000.0}]
}

proc timebetween_inhours {starttime endtime} {
	if {[catch {time_scan $endtime} endcode]} {return {}}
	if {[catch {time_scan $starttime} startcode]} {return {}}
	set diff [lmath_calc $endcode - $startcode]
	foreach {days miliseconds} $diff break
	expr {$days*24.0 + $miliseconds/3600000.0}
}

proc time_gt {time1 time2} {
	if {$time1 eq $time2} {return 0}
	if {[lsort -dict [list $time1 $time2]] eq [list $time1 $time2]} {return 0} else {return 1}
}

proc emptyifwrong {varVar} {
	upvar $varVar var
	if {![isdouble $var]} {set var ""}
	if {$var == 0} {set var ""}
	if {$var ne ""} {set var [format %.3f $var]}
}

proc time_seconds2duration {seconds} {
	set days [expr {int($seconds/86400)}]
	set miliseconds [expr {int(1000*($seconds - 86400*$days))}]
	set diff [list $days $miliseconds]
	timediff2duration $diff
}

proc time_seconds {diff} {
	foreach {days miliseconds} $diff break
	format %.3f [expr {$days*86400 + $miliseconds/1000.0}]
}

proc job_parse_log {job} {
	set submittime {} ; set starttime {} ; set endtime {} ; set duration {}
	set currentrun {} ; set currentsubmittime {}; set currentstarttime {} ; set currentjobid {} ; set currenthost {}
	set time_seconds ""
	set status unkown
	set logdata [split [file_read [job.file log $job]] \n]
	set failed 0
#	set tail [file tail $job]
	set tail .*
#	set poss [list_find -regexp $logdata submitted|running]
#	set logdata [lrange $logdata [lindex $poss end] end]
	set finishedlist {}
	set errorlist {}
	foreach line $logdata {
		# if {![regexp {submitted|starting|finished|failed|skipped|skipping} $line]} continue
		if {[regexp {^([0-9:. -]+)[ \t]-+ submitted .* \(run (.*)\) --} $line temp currentsubmittime currentrun]} {
			set status submitted
		} elseif {[regexp [subst -nocommands -nobackslashes {([0-9:. -]+)[ \t]starting ${tail} on (.*)($|:)}] $line temp currentstarttime currenthost]} {
			set status running
		} elseif {[regexp [subst -nocommands -nobackslashes {([0-9:. -]+)[ \t]starting ${tail}($|:)}] $line temp currentstarttime]} {
			set status running
		} elseif {[regexp [subst -nocommands -nobackslashes {([0-9:. -]+)[ \t]${tail} finished($|:)}] $line temp endtime]} {
			set run $currentrun
			set submittime $currentsubmittime
			set starttime $currentstarttime
			set status finished
			set host $currenthost
			set status finished
			lappend finishedlist [list $run $submittime $starttime $endtime $host]
		} elseif {
			[regexp [subst -nocommands -nobackslashes {([0-9:. -]+)[ \t]${tail} failed($|:)}] $line temp endtime] 
			|| [regexp [subst -nocommands -nobackslashes {([0-9:. -]+)[ \t]job ${tail} failed($|:)}] $line temp endtime]
		} {
			set run $currentrun
			set submittime $currentsubmittime
			set starttime $currentstarttime
			set status error
			set host $currenthost
			set status error
			lappend errorlist [list $run $submittime $starttime $endtime $host]
		} elseif {[regexp [subst -nocommands -nobackslashes {([0-9:. -]+)[ \t>-]+job ${tail} skipped($|:)}] $line temp skiptime]} {
			if {$status ne "finished"} {set status skipped}
		} elseif {[regexp [subst -nocommands -nobackslashes {([0-9:. -]+)[ \t]skipping ${tail}($|:)}] $line temp skiptime]} {
			if {$status ne "finished"} {set status skipped}
		}
	}
	if {$status eq "skipped"} {
		# if "skipped" it means that the result is there, so the job is finished, but
		# look for the last successful analysis to fill in times
		# if there is none, the log is incomplete (e.g. is finished outside of normal pipeline)
		# -> we should not check if still running or error though (as in next part)
		if {[llength $finishedlist]} {
			foreach {run submittime starttime endtime host} [lindex $finishedlist end] break
		} else {
			set submittime $skiptime
			set starttime $skiptime
			set endtime {}
			set host {}
		}
	} elseif {$status in "running submitted"} {
		set submittime $currentsubmittime
		set starttime $currentstarttime
		set host $currenthost
		set endtime {}
		set job_jid [job.file jid $job]
		if {![file exists $job_jid] || ![job_running [file_read $job_jid]]} {
			set status error
			if {[file exists [job.file err $job]]} {
				set endtime [clock format [file mtime [job.file err $job]] -format "%Y-%m-%d %H:%M:%S"]
				if {[string range $starttime 0 18] eq $endtime} {set endtime $starttime}
			} else {
				set endtime {}
			}
		} else {
			set status running
			set endtime {}
		}
	}
	if {![info exists run]} {set run $currentrun}
	if {$submittime eq ""} {set submittime $currentsubmittime}
	if {$starttime eq ""} {set submittime $currentstarttime}
	# putsvars submittime starttime endtime duration currentrun currentsubmittime currentstarttime
	if {$run eq ""} {
		set run [clock format [file mtime [job.file log $job]] -format "%Y-%m-%d_%H-%M-%S"]
	}
	set startcode [timescan starttime "error parsing [job.file log $job] (starttime)"]
	set endcode [timescan endtime "error parsing [job.file log $job] (endtime)"]
	if {$starttime ne ""} {
		if {$endtime ne ""} {
			set extratime {}
			set diff [lmath_calc $endcode - $startcode]
			set time_seconds [time_seconds $diff]
		} else {
			set endtime ""
			set extratime ...
			set endcode [time_scan [timestamp]]
			set diff [lmath_calc $endcode - $startcode]
		}
		set duration [timediff2duration $diff]$extratime
	}
	return [list $status $starttime $endtime $run $duration $submittime $time_seconds $host]
}

proc job_cleanmsg {msg} {
	string_change [string trim $msg] [list \t \\t \n \\n]
}

proc job_logfile_add {job jobid status {targets {}} {cores 1} {msg {}} {submittime {}} {starttime {}} {endtime {}}} {
	global cgjob
	if {[job_getinfo]} return
	if {![info exists cgjob(f_logfile)]} return
	set run [file tail [get cgjob(logfile) ""]]
	set cgjob(endtime) $endtime
	set msg [job_cleanmsg $msg]
	if {$starttime ne "" && $endtime ne ""} {
		set diff [lmath_calc [time_scan $endtime] - [time_scan $starttime]]
		set cgjob(totalduration) [lmath_calc $cgjob(totalduration) + $diff]
		set duration [timediff2duration $diff]
		set time_seconds [time_seconds $diff]
	} else {
		set duration ""
		set time_seconds ""
	}
	if {$cgjob(basedir) ne ""} {
		set pos [string length $cgjob(basedir)/]
		if {[string match $cgjob(basedir)/* $job]} {
			set job [string range $job $pos end]
		}
		set newtargets {}
		foreach target $targets {
			if {[string match $cgjob(basedir)/* $target]} {
				set target [string range $target $pos end]
			}
			lappend newtargets $target
		}
		set targets $newtargets
	}
	puts $cgjob(f_logfile) [join [list $job $jobid $status $submittime $starttime $endtime $duration $time_seconds $targets $msg $run $cores] \t]
	flush $cgjob(f_logfile)
	if {$status eq "error"} {set cgjob(status) error}
}

proc job_backup {file {rename 0}} {
	if {![job_file_or_link_exists $file]} return
	set num 1
	while 1 {
		if {![job_file_or_link_exists $file.old$num]} break
		incr num
	}
	if {$rename} {
		file rename -force -- $file $file.old$num
	} else {
		file copy $file $file.old$num
	}
}

proc job_to_old {file} {
	if {$::cgjob(dry)} {
		set cgjob_rm($file) old
	}
	if {![job_file_or_link_exists $file]} return
	file delete -force $file.old
	file rename -force -- $file $file.old
}

proc job_generate_code {job pwd adeps targetvars targets checkcompressed code} {
	set cmd ""
	set jobname [file tail $job]
	append cmd "file_add \{[job.file log $job]\} \"\[job_timestamp\]\\tstarting $jobname on \[exec hostname\]\"\n"
	append cmd "[list cd $pwd]\n"
	append cmd "[list set rootdir $pwd]\n"
	append cmd "[list set job $job]\n"
	append cmd "[list set jobname $jobname]\n"
	append cmd "[list set deps $adeps]\n"
	append cmd "[list set dep [lindex $adeps 0]]\n"
	append cmd "[list set checkcompressed $checkcompressed]\n"
	set num 1
	foreach dep $adeps {
		if {$num <= 10} {
			append cmd "[list set dep$num $dep]\n"
		}
		incr num
		if {$dep ne ""} {
			append cmd "if \{!\[[list job_file_or_link_exists $dep]\]\} \{error \"dependency $dep not found\"\}\n"
		}
	}
	set num 1
	foreach targetvar $targetvars {
		append cmd "[list set match$num $targetvar]\n"
		incr num
	}
	append cmd "[list set targets $targets]\n"
	append cmd "[list set target [lindex $targets 0]]\n"
	set num 1
	foreach target $targets {
		append cmd "[list set target$num $target]\n"
		incr num
		if {$num > 10} break
	}
	append cmd $code\n
	append cmd [string_change {
		set ok 1
		set errormsg {}
		cd {@PWD@}
		foreach target $targets {
			if {$checkcompressed} {
				set files [gzfiles $target]
			} else {
				set files [checkfiles $target]
			}
			if {[llength $files]} {
				file_add [job.file log $job] "[job_timestamp]\ttarget ok: $target"
			} else {
				set msg "target not found: $target"
				file_add [job.file log $job] "[job_timestamp]\t$msg"
				putslog $msg
				lappend errormsg $msg
				set ok 0
			}
		}
		if {$ok} {
			file_add [job.file log $job] "[job_timestamp]\t$jobname finished\n"
			catch {file delete [job.file pid $job]}
			catch {file delete [job.file jid $job]}
			set o [open [job.file err $job] a]
			puts $o "\nfinished [job_timestamp]"
			close $o
			catch {file rename -force -- [job.file err $job] [job.file ok $job]}
		} else {
			file_add [job.file log $job] "[job_timestamp]\tjob $jobname failed\n"
			error $errormsg
		}
	} [list @PWD@ $pwd]]
	return $cmd
}

# var targets contains list of all targets
# var target contains the first target
# var deps contains a list of all dependencies
# var dep contains the first element of the first dependency (so you do not have to do [lindex $deps 0] to get it)
# var match1, match2, ... contain the first, second, ... match (in parenthesis) in deps (foreach deps and plain deps)

proc job {jobname args} {
	global curjobid cgjob job_logdir_submit
	upvar job_logdir job_logdir
	if {![info exists job_logdir]} {
		error "The variable job_logdir is not set, This must be set before calling job, either by using the command job_logdir, or by setting the variable directly"
	}
	if {[info exists job_logdir_submit($jobname,$job_logdir)] && !$cgjob(resubmit)} {
		error "already submitted job $jobname for logdir $job_logdir"
	}
	set job_logdir_submit($jobname,$job_logdir) 1
	file mkdir $job_logdir
	if {![info exists cgjob(id)]} {set cgjob(id) 1}
	if {[llength $args] < 1} {error "wrong # args for target: must be job jobname -deps deps -targets targets -code code ..."}
	set rmtargets {}
	set pos 0
	set vars {}
	set procs {}
	set precode {}
	set skiplist {}
	set submitopts {}
	set checkcompressed 1
	set jobforce 0
	set optional 0
	set cores 1
	set len [llength $args]
	while {$pos < $len} {
		set key [lindex $args $pos]
		incr pos
		switch -- $key {
			-deps {
				set deps [lindex $args $pos]
				incr pos
			}
			-targets {
				set targets [lindex $args $pos]
				incr pos
			}
			-optional {
				set optional [lindex $args $pos]
				incr pos
			}
			-rmtargets {
				set rmtargets [lindex $args $pos]
				incr pos
			}
			-skip {
				lappend skiplist [lindex $args $pos]
				incr pos
			}
			-code {
				set code [lindex $args $pos]
				incr pos
			}
			-direct {
				set v [lindex $args $pos]
				if {$v ne "0"} {
					lappend submitopts -direct
				}
				if {[inlist {0 1} $v]} {incr pos}
			}
			-io {
				lappend submitopts -io [lindex $args $pos]
				incr pos
			}
			-cores {
				set cores [lindex $args $pos]
				lappend submitopts -cores $cores
				incr pos
			}
			-priority {
				set priority [lindex $args $pos]
				lappend submitopts -priority $priority
				incr pos
			}
			-mem {
				lappend submitopts -mem [lindex $args $pos]
				incr pos
			}
			-time {
				lappend submitopts -time [lindex $args $pos]
				incr pos
			}
			-hard {
				lappend submitopts -hard [lindex $args $pos]
				incr pos
			}
			-soft {
				lappend submitopts -soft [lindex $args $pos]
				incr pos
			}
			-vars {
				set vars [lindex $args $pos]
				incr pos
			}
			-procs {
				set procs [lindex $args $pos]
				incr pos
			}
			-precode {
				set precode [lindex $args $pos]
				incr pos
			}
			-checkcompressed {
				set checkcompressed [lindex $args $pos]
				incr pos
			}
			-force {
				set jobforce [lindex $args $pos]
				incr pos
			}
			-- break
			default {
				if {[string index $key 0] eq "-"} {
					error "unkown option $key for job, must be one of: -deps, -targets, -code, -vars, -procs, -rmtargets, -skip, -direct, -io, -cores, -precode, -checkcompressed"
				}
				break
			}
		}
	}
	if {$pos < $len} {
		set args [lrange $args [expr {$pos-1}] end]
		if {[llength $args] != 2 && [llength $args] != 3} {
			error "wrong # args for job: must be:\n job jobname -deps deps -targets targets -code code ... \nor\n job jobname options deps targets code"
		}
		lappend args {} {}
		set toget {}
		if {![info exists deps]} {lappend toget deps}
		if {![info exists targets]} {lappend toget targets}
		if {![info exists code]} {lappend toget code}
		foreach $toget $args break
	} else {
		if {![info exists deps]} {set deps {}}
		if {![info exists targets]} {set targets {}}
		if {![info exists code]} {set code {}}
	}
#	if {$targets eq ""} {
#		error "Each job must have targets (use -targets)"
#	}
	set level 1
	set edeps [job_expandvarslist $deps $level]
	set etargets [job_expandvarslist $targets $level]
	set ermtargets [job_expandvarslist $rmtargets $level]
	set eskip {}
	foreach skip $skiplist {
		lappend eskip [job_expandvarslist $skip $level]
	}
	set newcode {}
	if {[info exists ::defcompressionlevel]} {
		append newcode [list ::defcompressionlevel $::defcompressionlevel]\n
	}
	foreach var $vars {
		append newcode [list set $var [uplevel get $var]]\n
	}
	foreach proc $procs {
		append newcode [list proc $proc [info args $proc] [info body $proc]]\n
	}
	append newcode $code
	if {[get ::job_getinfo 0]} {
		# do not actually run if just gathering info
		job_process_getinfo $cgjob(id) $jobname $job_logdir [pwd] $edeps {} $etargets $eskip $checkcompressed $newcode $submitopts $ermtargets $precode $jobforce $optional $cores
		return
	}
	lappend cgjob(queue) [list $cgjob(id) $jobname $job_logdir [pwd] $edeps {} $etargets $eskip $checkcompressed $newcode $submitopts $ermtargets $precode $jobforce $optional $cores]
	incr cgjob(id)
	if {!$cgjob(debug)} {job_process}
}

proc job_init {args} {
	global cgjob cgjob_id cgjob_running job_logdir_submit cgjob_info cgjob_rm
	upvar job_logdir job_logdir
## job_init debugging code
#if {![info exists ::cgjob_initnum]} {
#	set ::cgjob_initnum 1
#} else {
#	incr ::cgjob_initnum
#}
#puts "----- job_init $::cgjob_initnum from proc [lindex [info level 1] 0]----"
#if {$::cgjob_initnum > 2} {
#	error "call job_init only once"
#}
	unset -nocomplain cgjob
	unset -nocomplain cgjob_id
	unset -nocomplain cgjob_rm
	unset -nocomplain cgjob_running
	unset -nocomplain cgjob_info
	unset -nocomplain job_logdir_submit
	set ::job_method_info {}
	set cgjob(hasargs) 0
	set cgjob(dry) 0
	set cgjob(debug) 0
	set cgjob(distribute) 0
	set cgjob(force) 0
	set cgjob(queue) {}
	set cgjob(id) 1
	set cgjob(resubmit) 1
	set cgjob(skipjoberrors) 0
	set cgjob(runcmd) [list $::genomecombdir/cg source -stack 1]
	set cgjob(logfile) {}
	set cgjob(starttime) {}
	set cgjob(endtime) {}
	set cgjob(totalduration) {0 0}
	set cgjob(cleanup) success
	set cgjob(removeold) 0
	set cgjob(cleanupshadowfiles) {}
	set cgjob(cleanupfiles) {}
	set cgjob(cleanupifemptyfiles) {}
	set cgjob(dmem) {}
	set cgjob(dtime) {}
	set job_logdir [file_absolute [pwd]/log_jobs]
	set cgjob(default_job_logdir) 1
	interp alias {} job_process {} job_process_direct
	interp alias {} job_runall {} job_runall_direct
	interp alias {} job_running {} job_running_direct
	interp alias {} job_wait {} job_wait_direct
	set ::job_getinfo 0
	job_args $args
}

proc job_curargs {} {
	global cgjob
	set temp ""
	foreach {opt field} {
		-distribute distribute -force force -dpriority priority -dqueue dqueue -dcleanup cleanup 
		-dmem {} -dtime {}
		-runcmd runcmd -skipjoberrors skipjoberrors
	} {
		if {[info exists cgjob($field)]} {
			lappend temp $opt $cgjob($field)
		}
	}
	return $temp
}

proc job_mempercore {mem threads} {
	if {[regexp {^([0-9]+)(.*)$} $mem temp memnum memunits]} {
		if {$memunits eq "G" || $memunits eq "g"} {
			set scale 1000.0
			set memunits M
		} else {
			set scale 1.0
		}
		set mem [expr {int(ceil($memnum * $scale /$threads))}]$memunits
	}
	return $mem
}

proc job_memgt {mem mem2} {
	if {$mem eq ""} {
		return 0
	}
	if {$mem2 eq ""} {
		return 1
	}
	set memnum $mem
	set memnum2 $mem2
	if {[regexp {^([0-9]+)(.*)$} $mem temp memnum memunits]} {
		if {$memunits eq "G" || $memunits eq "g"} {
			set memnum [expr {$memnum*1000.0}]
		}
	}
	if {[regexp {^([0-9]+)(.*)$} $mem2 temp2 memnum2 memunits2]} {
		if {$memunits2 eq "G" || $memunits2 eq "g"} {
			set memnum2 [expr {$memnum2*1000.0}]
		}
	}
	if {$memnum > $memnum2} {return 1} else {return 0}
}

proc job_cleanup_ifempty_add {args} {
	global cgjob
	foreach file $args {
		lappend cgjob(cleanupifemptyfiles) [file_absolute $file]
	}
}

proc job_cleanup_add {args} {
	global cgjob
	foreach file $args {
		lappend cgjob(cleanupfiles) [file_absolute $file]
	}
}

proc job_cleanup_add_shadow {args} {
	global cgjob
	foreach file $args {
		lappend cgjob(cleanupshadowfiles) [file_absolute $file]
	}
}

proc job_delete_ifempty {file {subdirs 0}} {
	if {[file isdir $file]} {
		set content [list_remove [glob -nocomplain $file/*] $file/shadow_source]
		if {$subdirs} {
			set empty 1
			foreach subdir $content {
				if {![file isdir $subdir] || [llength [glob -nocomplain $subdir/*]]} {
					set empty 0
				}
			}
			if {$empty} {set content {}}
		}
		if {![llength $content]} {
			shadow_delete $file
		}
	} else {
		shadow_delete $file
	}
	
}

proc job_cleanup {} {
	global cgjob
	foreach file $cgjob(cleanupshadowfiles) {
		catch {shadow_delete $file}
	}
	foreach file $cgjob(cleanupfiles) {
		catch {shadow_delete $file}
	}
	foreach file $cgjob(cleanupifemptyfiles) {
		job_delete_ifempty $file
	}
}

# maxsize: on linux filename size is usually limited to 255 characters.
# limit to 247 to allow for addition of the extension (.log, ...) and typedir
proc job_relfile2name {prefix file {maxsize 247}} {
	upvar job_logdir job_logdir
	if {[info exists job_logdir]} {
		set ffile [file_absolute $file]
		set dir $job_logdir
		while 1 {
			set dir [file dir $dir]
			if {$dir eq "/"} break
			set dirlength [string length $dir]
			incr dirlength -1
			if {[string range $ffile 0 $dirlength] eq $dir} {
				set file [string range $ffile [expr {$dirlength + 2}] end]
				if {$file eq ""} {set file [file tail $ffile]}
				break
			}
		}
	}
	set filepart [string_change $file {/ __ : _ \" _ \' _}]
	set len [string length $filepart]
	set start [expr {$len  - ($maxsize - [string length $prefix])}]
	return $prefix[string range $filepart $start end]
}

# This is used default in proc job to limit maxsize (for when job_relfile2name was not used)
proc shorten_filename {filename {maxsize 251}} {
	set tail [file tail $filename]
	set size [string length $tail]
	if {$size <= $maxsize} {return $filename}
	set dir [file dir $filename]
	if {$dir eq "."} {set dir {}} else {set dir $dir/}
	set pos [string first - $tail]
	set pos__ [string first __ $tail]
	if {$pos__ < $pos} {set pos [expr {$pos__ - 1}]}
	set prefix [string range $tail 0 $pos]
	set post [string range $tail [expr {$pos+1}] end]
	set size [string length $post]
	set start [expr {$size  - ($maxsize - [string length $prefix])}]
	return $dir$prefix[string range $post $start end]
}

proc clean_cmdline {args} {
	set cmdline "[list cd [pwd]] \; [list {*}$args]"
	regsub -all \n $cmdline " " cmdline
	return $cmdline
}

if {![info exists cgjob(distribute)]} {
	job_init
}
