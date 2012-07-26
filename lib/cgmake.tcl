# matches in the target can be used in the deps:
# if target is dir/test_2, that matches target pattern (.*)/test_(.*), then the following variables can be used in deps:
# target = dir/test_2
# target1 = dir
# target2 = 2
#
# dependencies between braces () are optional (braces must be at start and end of dependency)
# 
# variables in vars or used in deps are automatically made (targets) using direct
# if the variable is preceeded by a _ (e.g. $_var), the dependency will be expanded:
# var will be used as target, and is expected to contain a list of values
# the dependencies will be expanded: a list with a dependency for each element of the var list
proc target {args} {
	global cgmakedata
	if {[llength $args] < 3} {error "wrong # args for target: must be target targetname options"}
	set targetname [lindex $args 0]
	set target [lindex $args 1]
	set pos 2
	set deps {}
	set vars {}
	set precode {}
	set code {}
	set submitopts {}
	while {$pos < [llength $args]} {
		set key [lindex $args $pos]
		incr pos
		switch -- $key {
			-deps {
				set deps [lindex $args $pos]
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
			-- break
			default {
				if {[string index $key 0] eq "-"} {
					error "unkown option $key for target, must be one of: -deps, -direct, -io"
				} else {
					error "wrong # args for target: must be target targetname options"
				}
				break
			}
		}
	}
	lappend cgmakedata($target) [list $targetname $target $deps $vars $precode $code $submitopts]
}

proc cgmake_args {args} {
	global cgmakeargs
	if {![info exists cgmakeargs(distribute)]} {
		set cgmakeargs(distribute) 0
	}
	if {![llength $args]} {return {}}
	set newargs {}
	set pos 0
	while {$pos < [llength $args]} {
		set key [lindex $args $pos]
		incr pos
		switch -- $key {
			-d - --distribute {
				set cgmakeargs(distribute) [lindex $args $pos]
				incr pos
			}
			-- break
			default {
				lappend newargs $key
			}
		}
	}
	lappend newargs {*}[lrange $args $pos end]
}

proc cgmake_replacetarget {string target targetvars} {
	set poss [regexp -all -inline -indices {\$\{?target[0-9]*\}?} $string]
	list_foreach {s e} [list_reverse $poss] {
		set t [string range $string $s $e]
		if {[regexp {[0-9]+} $t num]} {
			incr num -1
			set string [string_replace $string $s $e [lindex $targetvars $num]]
		} else {
			set string [string_replace $string $s $e $target]
		}
	}
	return $string
}

proc cgmake_expandvars {string} {
	set code "proc cgmake_testvars \{string\} \{"
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
			eval "$code\nsubst \$string\n\}"
			if {![catch {cgmake_testvars $string} result]} {
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
			set status [cgmaketarget $var newids 1]
			if {[regexp {^cannot make} $status]} {
				putslog $status
				error $status
			}
			if {![info exists ::$var]} {
				error "cannot make $var: make finished, but var is not defined"
			}
			if {$expand} {
				set valuelist [get ::$var]
				foreach value $valuelist {
					lappend resultlist "$code\n[list set _$var $value]" 0 {}
				}
			} else {
				append code "\n[list set $var [get ::$var]]"
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

proc cgmake_expandvarslist {list} {
	set result {}
	foreach string $list {
		lappend result {*}[cgmake_expandvars $string]
	}
	return $result
}

proc cgmakewait {ids} {
}

proc putslog {args} {
	global loglevel
	foreach message $args {
		puts stderr "[get loglevel ""]$message"
	}
}

proc logleveldown {} {
	append ::loglevel "    "
}

proc loglevelup {} {
	set ::loglevel [string range $::loglevel 0 end-4]
}

proc cgmaketarget {target newidsVar {direct 0}} {
	# for later if direct=1, run direct
	global cgmakedata cgmakeargs cgmakeids cgmakeroot
	upvar $newidsVar newids
	putslog "================================================================="
	putslog "Target $target"
	cd $cgmakeroot
	set newids {}
	set target [uplevel #0 [list subst -nobackslashes -nocommands $target]]
	if {[llength [gzfiles $target]]} {
		putslog "file $target exists"
		return ok
	} elseif {[info exists cgmakeids($target)]} {
		if {$cgmakeids($target) eq ""} {
			putslog "target $target already done"
		} else {
			putslog "target $target already submitted (id $cgmakeids($target))"
			lappend newids $cgmakeids($target)
		}
		return ok
	} elseif {[info exists ::$target]} {
		putslog "variable $target exists"
		return ok
	}
	logleveldown
	putslog "Finding rules for target $target"
	set checklist {}
	if {[info exists cgmakedata($target)]} {
		lappend checklist {*}$cgmakedata($target)
	}
	foreach pattern [array names cgmakedata] {
		set spattern [uplevel #0 [list subst -nobackslashes -nocommands $pattern]]
		if {[regexp ^$spattern\$ $target]} {
			lappend checklist {*}$cgmakedata($pattern)
		}
	}
	if {![llength $checklist]} {
		putslog "cannot make $target: no rule found"
		loglevelup
		return "cannot make $target: no rule found"
	}
	list_foreach {targetname pattern deps vars precode code submitopts} $checklist {
		putslog "Test rule $targetname"
		putslog "     pattern: $pattern"
		putslog "     dep: [join $deps "\n    dep: "]"
		set targetvars {}
		set pattern [uplevel #0 [list subst -nobackslashes -nocommands $pattern]]
		if {$target ne $pattern} {
			set targetvars [lrange [regexp -all -inline ^$pattern\$ $target] 1 end]
			if {[llength $targetvars]} {
				set num 1
				foreach v $targetvars {
					set ::target$num  $v
					incr num
				}
			}
		}
		set ids {}
		set depsok 1
		set deps [cgmake_replacetarget $deps $target $targetvars]
		set error [catch {cgmake_expandvarslist $deps} deps]
		if {$error} {
			if {![regexp {^cannot make} $deps]} {error $deps $::errorInfo}
			# one of the deps cannot be made, try next rule
			set depsok 0
			continue
		}
		set vars [cgmake_replacetarget $vars $target $targetvars]
		set error [catch {cgmake_expandvarslist $vars} vars]
		if {$error} {
			if {![regexp {^cannot make} $deps]} {error $deps $::errorInfo}
			# one of the deps cannot be made, try next rule
			set depsok 0
			continue
		}
		if {$precode ne ""} {
			putslog "Running precode: $precode"
			set precode [cgmakecode $precode $target $targetvars $deps $vars]
			proc cgmake_temp {} $precode
			uplevel #0 cgmake_temp
		}
		foreach dep $deps {
			if {[string index $dep 0] eq "\(" && [string index $dep end] eq "\)"} {
				set opt 1
				set dep [string range $dep 1 end-1]
			} else {
				set opt 0
			}
			set status [cgmaketarget $dep newids]
			if {!$opt && [regexp {^cannot make} $status]} {
				set depsok 0
				break
			}
			if {[llength $newids]} {
				lappend ids {*}$newids
			}
		}
		if {$depsok} break
		set depsok 0
	}
	if {!$depsok} {
		putslog "cannot make $target: missing dependencies"
		loglevelup
		return "cannot make $target: missing dependencies"
	}
	loglevelup
	putslog "making $target (rule $targetname)"
	set ecode [cgmakecode $code $target $targetvars $deps $vars]
	if {$cgmakeargs(distribute) == 0} {
#		incr cgmakeids()
#		set cgmakeids($target) $cgmakeids()
		proc cgmake_temp {} $ecode
		uplevel #0 cgmake_temp
		if {[llength [gzfiles $target]]} {
			putslog "file $target ok"
			return ok
		} elseif {[info exists ::$target]} {
			putslog "variable $target ok"
			return ok
		} else {
			putslog "cannot make $target: ERROR rule $targetname did not actually make target"
			return "cannot make $target: ERROR rule $targetname did not actually make target"
		}
#		set cgmakeids($target) {}
		putslog "Made $targetname ($target)\n"
	} else {
		error "other than -d 0 not implemented yet"
	}
	return ok
}

proc cgmakecode {code target targetvars deps vars} {
	global cgmakeroot
	set ecode {}
	append ecode [list cd $cgmakeroot]\n
	append ecode [list set deps $deps]\n
	append ecode [list set dep [lindex $deps 0]]\n
	append ecode [list set target $target]\n
	set num 1
	foreach v $targetvars {
		append ecode [list set target$num $v]\n
		incr num
	}
	foreach var $vars {
		append ecode [list set $var [get ::$var]]\n
	}
	append ecode $code
}

proc cgmake_clear {} {
	global cgmakedata cgmakeargs cgmakeids
	unset -nocomplain cgmakedata
	unset -nocomplain cgmakeargs
	unset -nocomplain cgmakeids
}

proc cgmake_targets {} {
	global cgmakedata
	array names cgmakedata
}

proc cgmake {args} {
	global cgmakeargs cgmakeids cgmakeroot loglevel
	set loglevel ""
	unset -nocomplain cgmakeids
	set cgmakeroot [pwd]
	set cgmakeids() 0
	set args [cgmake_args {*}$args]
	if {![llength $args]} {set args all}
	set targetsnotok {}
	foreach target $args {
		if {[string index $target 0] eq "\(" && [string index $target end] eq "\)"} {
			set opt 1
			set dep [string range $target 1 end-1]
		} else {
			set opt 0
		}
		set status [cgmaketarget $target newids]
		if {!$opt && [regexp {^cannot make} $status]} {
			lappend targetsnotok $target
		}
	}
	if {![llength $targetsnotok]} {
		putslog "cgmake finished"
	} else {
		putslog "cgmake unfinished, missed targets:\n  [join $targetsnotok "\n  "]"
	}
}

if 0 {

	cgmake_clear
	unset -nocomplain names

	set cgmakefile ../cgmake.tcl ../tests/data test

	mkdir ~/dev/genomecomb/tests/tmp
	cd ~/dev/genomecomb/tests/tmp
	rm ~/dev/genomecomb/tests/tmp/*
	cp cgmaketest data/project.tsv ~/dev/genomecomb/tests/tmp/
	cg make -cgmake ../cgmaketest -projectfile ../data/project.tsv test

}

proc cg_make {args} {
	cgmake_clear
	
}
