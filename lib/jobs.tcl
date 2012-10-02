proc regexp2glob {pattern} {
	regsub -all {[\[\{(][^\]\})]*[\]\})]} $pattern {*} glob
	regsub -all {\\.} $glob {*} glob
	regsub -all {[*+?.]+} $glob {*} glob
	return $glob
}

proc job_args {jobargs} {
	global cgjob
	if {![info exists cgjob(distribute)]} {
		set cgjob(distribute) 0
	}
	if {![info exists cgjob(force)]} {
		set cgjob(force) 0
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
				set cgjob(distribute) [lindex $jobargs $pos]
				incr pos
			}
			-- break
			default {
				lappend newargs $key
			}
		}
	}
	lappend newargs {*}[lrange $jobargs $pos end]
}

proc job_init {args} {
	global cgjob job_logdir
	set cgjob(distribute) 0
	set cgjob(force) 0
	set cgjob(queue) {}
	set cgjob(id) 1
	job_logdir [file normalize [pwd]/log_jobs]
}

job_init

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
			set cmd "$code\nsubst -nobackslashes \$string\n\}"
			eval $cmd
			if {![catch {job_testvars $string} result]} {
				lappend resultlist $code 1 $result
				continue
			}
			if {![regexp {can't read "(.*)": no such variable} $result temp var]} {
				putslog "cannot make $var: ERROR: $result"
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

proc job_finddep {pattern} {
	return [lsort -dict [gzfiles $pattern]]
}

proc job_findregexpdep {pattern} {
	global job_data job_files
	set glob [regexp2glob $pattern]
	set files {}
	foreach file [lsort -dict [gzfiles $glob]] {
		if {[regexp ^$pattern\$ $file]} {lappend files $file}
	}
	return $files
}

# set string {abc\1def\\1ghi\\\1jkl\\\\1s}
proc job_targetreplace {string targetvars} {
	if {![llength $targetvars]} {return $string}
	list_foreach {b e} [list_reverse [regexp -inline -all -indices {\\+[0-9]+} $string]] {
		regexp {(\\+)([0-9]+)} [string range $string $b $e] temp pre num
		set len [string length $pre]
		if {![expr {$len%2}]} continue
		incr num -1
		set string [string_replace $string [expr {$b+$len-1}] $e [lindex $targetvars $num]]
	}
	return $string
}

proc job_targetsreplace {list targetvars} {
	set result {}
	foreach string $list {
		lappend result [job_targetreplace $string $targetvars]
	}
	return $result
}

proc job_finddeps {deps targetvarsVar {ftargetvars {}}} {
	upvar $targetvarsVar targetvars
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
			set files [job_findregexpdep $pattern]
		} else {
			set files [job_finddep $pattern]
		}
		if {![llength $files]} {
			if {!$opt} {
				error "missing dependency $pattern"
			}
			continue
		}
		lappend finaldeps {*}$files
		foreach file $files {
			lappend targetvars {*}[lrange [regexp -all -inline ^$pattern\$ $file] 1 end]
		}
	}
	return $finaldeps
}

# dependencies between braces () are optional (braces must be at start and end of dependency)

proc job_checktarget {target {newidsVar {}}} {
	global cgjob_ids
	if {$newidsVar ne ""} {
		upvar $newidsVar newids
	}
	if {[llength [gzfiles $target]]} {
		job_log "target $target ok"
		return 1
	} elseif {[info exists cgjob_ids($target)]} {
		if {$cgjob_ids($target) eq ""} {
			job_log "target $target already done"
			return 1
		} else {
			job_log "target $target already submitted (id $cgjob_ids($target))"
			set newids $cgjob_ids($target)
			return 2
		}
	} elseif {[info exists ::$target]} {
		job_log "variable $target exists (= [get ::$target])"
		return 1
	}
	return 0
}

proc job_checktargets {targets {runningVar {}}} {
	if {$runningVar ne ""} {
		upvar $runningVar running
		set running {}
	}
	set ok 1
	foreach target $targets {
		set check [job_checktarget $target]
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
	global cgjob_ids
	set targets {}
	set ok 1
	foreach ptarget $ptargets {
		lappend targets {*}[job_finddep $ptarget]
	}
	return $targets
}

proc job_logdir {logdir} {
	global job_log
	if {![info exists job_log(file)]} {set job_log(file) joblog.txt}
	set job_log(file) [file normalize $logdir/[file tail $job_log(file)]]
	file mkdir $logdir
	catch {close $job_log(f)}
	set job_log(f) [open $job_log(file) a+]
}

proc job_logname {name} {
	global job_log
	set name log_[string_change $name {/ __ : _ \" _ \' _}].txt
	set job_log(file) [file normalize [file dir $job_log(file)]/$name]
	catch {close $job_log(f)}
	set job_log(f) [open $job_log(file) a+]
}

proc job_log {args} {
	global job_log
	foreach message $args {
		puts $job_log(f) "[timestamp] $message"
		puts stderr "[timestamp] $message"
	}
	flush $job_log(f)
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

proc job_process {} {
	global cgjob job_deps 
	set jobroot [pwd]
	foreach line $cgjob(queue) {
putsvars line
		foreach {jobid mjobname pwd deps foreach ftargets fptargets fskip code submitopts} $line break
		cd $pwd
		set newids {}
		job_logname $mjobname
		# check foreach deps, skip if not fullfilled
		if {[llength $foreach]} {
			if {[catch {job_finddeps $foreach ftargetvars} fadeps]} {
				job_log "error in foreach dependencies for $mjobname: $fadeps"
				job_log "job $mjobname failed"
				continue
			}
			set fadeps [list_concat $fadeps]
			set ftargetvars [list_concat $ftargetvars]
		} else {
			set fadeps {{}}
			set ftargetvars {{}}
		}
		foreach fdep $fadeps ftargetvar $ftargetvars {
			cd $pwd
			if {$fdep eq ""} {
				set jobname $mjobname
			} else {
				set jobname $mjobname-$fdep
			}
			job_logname $jobname
			job_log "==================== $jobname ===================="
			# check deps, skip if not fullfilled
			if {[catch {job_finddeps $deps newtargetvars $ftargetvar} newadeps]} {
				job_log "error in dependencies for $jobname: $newadeps"
				job_log "job $jobname failed"
				continue
			}
			set targetvars $ftargetvar
			lappend targetvars {*}$newtargetvars
			set adeps $fdep
			lappend adeps {*}$newadeps
			# check targets, if already done or running, skip
			set run 0
			if {!$cgjob(force) && [llength $fskip]} {
				set skip [job_targetsreplace $fskip $targetvars]
				if {[llength $skip] && [job_checktargets $skip running]} {
					job_log "skipping $jobname: skip targets already completed or running"
					continue
				}
			}
			set targets [job_targetsreplace $ftargets $targetvars]
			if {![job_checktargets $targets running]} {
				set run 1
			}
			set ptargets [job_targetsreplace $fptargets $targetvars]
			if {[llength $ptargets] && ![llength [job_findptargets $ptargets]]} {
				set run 1
			}
			if {$cgjob(force)} {
				if {[llength $running]} {
					error "cannot force job with still running tasks ([join $running ,])"
				}
				foreach target [list_concat $targets $ptargets] {
					job_backup $target 1
				}
				set run 1
			}
			if {!$run} {
				job_log "skipping $jobname: targets already completed or running"
				continue
			}
			# run code
			set cmd "proc job_run {} \{\n"
			append cmd "[list cd $pwd]\n"
			append cmd "[list set deps $adeps]\n"
			append cmd "[list set dep [lindex $adeps 0]]\n"
			set num 1
			foreach targetvar $targetvars {
				append cmd "[list set dep$num $targetvar]\n"
				incr num
			}
			append cmd "[list set targets $targets]\n"
			append cmd "[list set target [lindex $targets 0]]\n"
			append cmd $code\n\}
			set ok 1
			if {[catch {eval $cmd} result]} {
				set ok 0
				job_log "error creating $jobname: $result"
			}
			if {[catch {job_run} result]} {
				set ok 0
				job_log "error running $jobname: $result"
			}
			# check if targets are ok
			if {![job_checktargets $targets]} {
				set ok 0
				job_log "job $jobname failed: missing targets"
			}
			if {[llength $ptargets] && ![llength [job_findptargets $ptargets]]} {
				set ok 0
				job_log "job $jobname failed: missing ptargets"
			}
			if {$ok} {
				job_log "job $jobname success"
			}
		}
	}
	cd $jobroot
	set cgjob(queue) {}
}

# var targets contains list of all targets
# var target contains the first target
# var deps contains a list of all dependencies
# var dep contains the first element of the first dependency (so you do not have to do [lindex $deps 0] to get it)

proc job {jobname args} {
	global curjobid cgjob
	if {![info exists cgjob(id)]} {set cgjob(id) 1}
	if {[llength $args] < 1} {error "wrong # args for target: must be job jobname -deps deps -targets targets -code code ..."}
	set pos 0
	set foreach {}
	set deps {}
	set vars {}
	set precode {}
	set code {}
	set targets {}
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
				} else {
					error "wrong # args for target: must be job submit jobname -deps deps -targets targets -code code ..."
				}
				break
			}
		}
	}
	if {$pos < $len} {
		error "wrong # args for target: must be job submit jobname -deps deps -targets targets -code code ..."
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
	lappend cgjob(queue) [list $cgjob(id) $jobname [pwd] $edeps $eforeach $etargets $eptargets $eskip $newcode $submitopts $precode]	
	incr cgjob(id)
	job_process
}

