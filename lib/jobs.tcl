
# set pattern (a|(bc|de))X(c|d)
# set pattern adfg(ab)
proc regexp2glob {pattern} {
	set glob $pattern
	while {[regsub -all {\([^\(\)]*[|][^\(\)]*\)} $glob {*} glob]} {}
	regsub -all {\(\[^\(\)]*|[^\(\)]*\)} $pattern {*} glob
	regsub -all {\[[^]]*\]} $glob {*} glob
	regsub -all {\{[^\}]*\}} $glob {*} glob
	regsub -all {\\.} $glob {*} glob
	regsub -all {[*+?.()]+} $glob {*} glob
	return $glob
}

proc job_distribute {type} {
	global cgjob
	set cgjob(distribute) $type
	if {[isint $cgjob(distribute)]} {
		if {$cgjob(distribute) <= 1} {
			set target direct
		} else {
			set target distr
		}
	} else {
		set target $type
	}
	auto_load job_process_$target
	interp alias {} job_process {} job_process_$target
	interp alias {} job_wait {} job_process_${target}_wait
}

proc job_args {jobargs} {
	global cgjob
	if {![info exists cgjob(distribute)]} {
		set cgjob(distribute) 0
	}
	if {![info exists cgjob(force)]} {
		set cgjob(force) 0
	}
	if {![info exists cgjob(silent)]} {
		set cgjob(silent) 0
	}
	if {![info exists cgjob(debug)]} {
		set cgjob(debug) 0
	}
	if {![info exists cgjob(resubmit)]} {
		set cgjob(resubmit) 0
	}
	if {![llength $jobargs]} {return {}}
	set newargs {}
	set pos 0
	while {$pos < [llength $jobargs]} {
		set key [lindex $jobargs $pos]
		incr pos
		switch -- $key {
			-f - -force - --force {
				set cgjob(force) [lindex $jobargs $pos]
				incr pos
			}
			-d - -distribute - --distribute {
				job_distribute [lindex $jobargs $pos]
				incr pos
			}
			-s - -silent - --silent {
				set cgjob(silent) 1
			}
			-debug {
				set cgjob(debug) 1
			}
			-resubmit {
				set cgjob(resubmit) 1
			}
			-- break
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
			upvar $level $var value
			if {![info exists value]} {
				error "Error: can't read \"$var\": no such variable $var"
			}
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

proc job_expandvarslist {list {level 1}} {
	set result {}
	incr level
	set list [string_change $list [list \\ \\\\]]
	foreach string $list {
		lappend result {*}[job_expandvars $string $level]
	}
	return $result
}

proc job_finddep {pattern idsVar timeVar} {
	global cgjob_id cgjob_ptargets
	upvar $idsVar ids
	upvar $timeVar time
	set pattern [file normalize $pattern]
	set ptargethits [array names cgjob_ptargets $pattern]
	if {[llength $ptargethits]} {
		error "ptargets hit $pattern: wait till ptarget deps have finished"
	}
	set files [lsort -dict [gzfiles $pattern]]
	foreach file $files {
		maxfiletime $file time
	}
	lappend ids {*}[list_fill [llength $files] {}]
	foreach file [array names cgjob_id $pattern] {
		if {[inlist $files $file]} {
			unset cgjob_id($file)
			continue
		}
		lappend files $file
		lappend ids $cgjob_id($file)
		set time now
	}
	return $files
}

proc maxfiletime {file timeVar} {
	upvar $timeVar time
	if {$time eq "now"} {return $time}
	set ftime [file mtime $file]
	if {$ftime > $time} {set time $ftime}
}

proc job_findregexpdep {pattern idsVar timeVar} {
	global cgjob_id cgjob_ptargets
	upvar $idsVar ids
	upvar $timeVar time
	set pattern [file normalize $pattern]
	set glob [regexp2glob $pattern]
	if {[llength [array names cgjob_ptargets $glob]]} {
		error "ptargets hit $pattern: wait till ptarget deps have finished"
	}
	set files {}
	foreach file [lsort -dict [gzfiles $glob]] {
		if {[regexp ^$pattern\$ $file]} {
			maxfiletime $file time
			lappend files $file
			lappend ids {}
		}
	}
	foreach file [array names cgjob_id [file normalize $glob]] {
		if {![regexp ^[file normalize $pattern]\$ $file]} continue
		if {[inlist $files $file]} {
			unset cgjob_id($file)
			continue
		}
		lappend files $file
		lappend ids $cgjob_id($file)
		set time now
	}
	return $files
}

# dependencies between braces () are optional (braces must be at start and end of dependency)
proc job_finddeps {job deps targetvarsVar targetvarslist idsVar timeVar {ftargetvars {}}} {
	upvar $idsVar ids
	if {$targetvarsVar ne ""} {
		upvar $targetvarsVar targetvars
	}
	upvar $timeVar time
	set time 0
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
		if {[string index $pattern 0] eq "^" && [string index $pattern end] eq "\$"} {
			set pattern [string range $pattern 1 end-1]
			set files [job_findregexpdep $pattern ids time]
		} else {
			set files [job_finddep $pattern ids time]
		}
		if {![llength $files]} {
			if {!$opt} {
				error "missing dependency $pattern"
			} else {
				job_lognf $job "missing optional dependency $pattern"
				lappend finaldeps {}
			}
			continue
		} else {
			job_lognf $job "dependency ok ($pattern): $files"
		}
		lappend finaldeps {*}$files
		foreach file $files {
			set targets [lrange [regexp -all -inline ^[file normalize $pattern]\$ $file] 1 end]
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

# set string {abc\1def\\1ghi\\\1jkl\\\\1s}
proc job_targetreplace {string targetvars} {
	set targetvarslen [llength $targetvars]
	list_foreach {b e} [list_reverse [regexp -inline -all -indices {\\+[0-9]+} $string]] {
		regexp {(\\+)([0-9]+)} [string range $string $b $e] temp pre num
		set len [string length $pre]
		if {![expr {$len%2}]} continue
		incr num -1
		if {!$targetvarslen} {set value *} else {set value [lindex $targetvars $num]}
		set string [string_replace $string [expr {$b+$len-1}] $e $value]
	}
	return $string
}

proc job_targetsreplace {list targetvars} {
	set result {}
	foreach string $list {
		lappend result [file normalize [job_targetreplace $string $targetvars]]
	}
	return $result
}

proc job_checktarget {job target time {newidsVar {}}} {
	global cgjob_id
	if {$newidsVar ne ""} {
		upvar $newidsVar newids
	}
	set files [gzfiles $target]
	if {[llength $files]} {
		if {$time ne "now"} {
			foreach file $files {
				if {[file mtime $file] < $time} {
					set time now
					break
				}
			}
		}
		if {$time eq "now"} {
			job_lognf $job "target older than dep (removing): $target"
			# foreach file $files {
			# 	job_backup $file 1
			# }
			file delete {*}$files
			return 0
		} else {		
			job_lognf $job "target ok: $target"
			unset -nocomplain cgjob_id($target)
			return 1
		}
	} elseif {[info exists cgjob_id($target)] && $cgjob_id($target) != "q"} {
		job_lognf $job "target already submitted/running (id $cgjob_id($target)): $target "
		set newids $cgjob_id($target)
		return 2
	}
	job_lognf $job "target missing: $target"
	return 0
}

proc job_checktargets {job targets time {runningVar {}}} {
	if {$runningVar ne ""} {
		upvar $runningVar running
		set running {}
	}
	set ok 1
	foreach target $targets {
		set check [job_checktarget $job $target $time]
		if {!$check} {
			set ok 0
		}
		if {$check == 2} {
			lappend running $target
		}
	}
	return $ok
}

proc job_findptargets {ptargets} {
	global cgjob_id
	set targets {}
	set ok 1
	foreach pattern $ptargets {
		if {[string index $pattern 0] eq "^" && [string index $pattern end] eq "\$"} {
			set pattern [string range $pattern 1 end-1]
			set glob [regexp2glob $pattern]
			set files {}
			foreach file [lsort -dict [gzfiles $glob]] {
				if {[regexp ^$pattern\$ $file]} {
					lappend files $file
					lappend ids {}
				}
			}
		} else {
			set files [lsort -dict [gzfiles $pattern]]
		}
		lappend targets {*}$files
	}
	return $targets
}

proc job_logdir {{logdir {}}} {
	upvar job_logdir job_logdir
	if {$logdir eq ""} {
		set logdir [file join [pwd] log_jobs]
	}
	set job_logdir [file normalize $logdir]
}

proc job_logname {job_logdir name} {
	set name [string_change $name {/ __ : _ \" _ \' _}]
	return [file normalize $job_logdir/$name]
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
		set cgjob(f,$job) [open $job.log a]
	}
	set f $cgjob(f,$job)
	set log [join $cgjob(buffer,$job) \n]
	puts $f $log
	if {!$cgjob(silent)} {puts stderr $log}
	set cgjob(buffer,$job) {}
	flush $f
}

proc job_lognf {job args} {
	global cgjob
	foreach message $args {
		lappend cgjob(buffer,$job) "[job_timestamp]\t$message"
	}
}

proc job_logclose {job args} {
	global cgjob
	if {![info exists cgjob(f,$job)]} {
		set cgjob(f,$job) [open $job.log a]
	}
	job_log $job
	set f $cgjob(f,$job)
	catch {close $f}
	unset cgjob(f,$job)
	unset cgjob(buffer,$job)
}

proc job_backup {file {rename 0}} {
	if {![file exists $file]} return
	set num 1
	while 1 {
		if {![file exists $file.old$num]} break
		incr num
	}
	if {$rename} {
		file rename $file $file.old$num
	} else {
		file copy $file $file.old$num
	}
}

proc job_generate_code {job pwd adeps targetvars targets ptargets code} {
	set cmd ""
	set jobname [file tail $job]
	append cmd "file_add \{$job.log\} \"[job_timestamp]\\tstarting $jobname\"\n"
	append cmd "[list cd $pwd]\n"
	append cmd "[list set rootdir $pwd]\n"
	append cmd "[list set job $job]\n"
	append cmd "[list set jobname $jobname]\n"
	append cmd "[list set ptargets $ptargets]\n"
	append cmd "[list set deps $adeps]\n"
	append cmd "[list set dep [lindex $adeps 0]]\n"
	set num 1
	foreach dep $adeps {
		append cmd "[list set dep$num $dep]\n"
		incr num
		if {$num > 10} break
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
		cd {@PWD@}
		foreach target @TARGETS@ {
			set files [gzfiles $target]
			if {[llength $files]} {
				file_add {@JOB@.log} "[job_timestamp]\ttarget ok: $target"
			} else {
				file_add {@JOB@.log} "[job_timestamp]\ttarget not found: $target"
				set ok 0
			}
		}
		if {[llength @PTARGETS@]} {
			if {[llength [job_findptargets @PTARGETS@]]} {
				file_add {@JOB@.log} "[job_timestamp]\tptargets ok"
			} else {
				set ok 0
				file_add {@JOB@.log} "[job_timestamp]\tmissing ptargets"
			}
		}
		if {$ok} {
			file_add {@JOB@.log} "[job_timestamp]\tfinished @JOBNAME@\n"
		} else {
			file_add {@JOB@.log} "[job_timestamp]\tfailed @JOBNAME@\n"
		}
	} [list @PWD@ $pwd @JOB@ $job @JOBNAME@ $jobname @TARGETS@ [list $targets] @PTARGETS@ [list $ptargets]]]
	append cmd "file_write $job.finished \[timestamp\]\n"
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
	if {![info exists cgjob(id)]} {set cgjob(id) 1}
	if {[llength $args] < 1} {error "wrong # args for target: must be job jobname -deps deps -targets targets -code code ..."}
	set pos 0
	set foreach {}
	set vars {}
	set precode {}
	set skip {}
	set ptargets {}
	set submitopts {}
	set len [llength $args]
	while {$pos < $len} {
		set key [lindex $args $pos]
		incr pos
		switch -- $key {
			-deps {
				set deps [lindex $args $pos]
				incr pos
			}
			-foreach {
				set foreach [lindex $args $pos]
				incr pos
			}
			-targets {
				set targets [lindex $args $pos]
				incr pos
			}
			-skip {
				set skip [lindex $args $pos]
				incr pos
			}
			-ptargets {
				set ptargets [lindex $args $pos]
				incr pos
			}
			-code {
				set code [lindex $args $pos]
				incr pos
			}
			-direct {
				lappend submitopts -direct
			}
			-io {
				lappend submitopts -io [lindex $args $pos]
				incr pos
			}
			-vars {
				set vars [lindex $args $pos]
				incr pos
			}
			-precode {
				set precode [lindex $args $pos]
				incr pos
			}
			-- break
			default {
				if {[string index $key 0] eq "-"} {
					error "unkown option $key for target, must be one of: -deps, -targets, -code, -direct, -io"
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
	if {$targets eq ""} {
		error "Each job must have targets (use -targets)"
	}
	set edeps [job_expandvarslist $deps 1]
	set eforeach [job_expandvarslist $foreach 1]
	set etargets [job_expandvarslist $targets 1]
	set eskip [job_expandvarslist $skip 1]
	set eptargets [job_expandvarslist $ptargets 1]
	set newcode {}
	foreach var $vars {
		append newcode [list set $var [uplevel get $var]]\n
	}
	append newcode $code
	lappend cgjob(queue) [list $cgjob(id) $jobname $job_logdir [pwd] $edeps $eforeach {} $etargets $eptargets $eskip $newcode $submitopts $precode]	
	incr cgjob(id)
	if {!$cgjob(debug)} job_process
}

proc job_init {args} {
	global cgjob cgjob_id cgjob_running cgjob_ptargets job_logdir_submit
	upvar job_logdir job_logdir
	unset -nocomplain cgjob
	unset -nocomplain cgjob_id
	unset -nocomplain cgjob_running
	unset -nocomplain cgjob_ptargets
	unset -nocomplain job_logdir_submit
	set cgjob(debug) 0
	set cgjob(distribute) 0
	set cgjob(force) 0
	set cgjob(queue) {}
	set cgjob(id) 1
	set cgjob(resubmit) 0
	set job_logdir [file normalize [pwd]/log_jobs]
	interp alias {} job_process {} job_process_direct
	interp alias {} job_wait {} job_process_direct_wait
	job_args $args
}

job_init

