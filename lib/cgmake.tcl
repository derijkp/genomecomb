proc target {args} {
	global cgmakedata
	if {[llength $args] < 2} {exiterror "wrong # args for target: must be target targetname options"}
	set target [lindex $args 0]
	set pos 1
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
					exiterror "unkown option $key for target, must be one of: -deps, -direct, -io"
				} else {
					exiterror "wrong # args for target: must be target targetname options"
				}
				break
			}
		}
	}
	lappend cgmakedata($target) [list $target $deps $vars $precode $code $submitopts]
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

proc cgmake_expanddeps {deps {targetvars {}}} {
	set newdeps {}
	set result {}
	foreach dep $deps {
		if {[llength $targetvars]} {
			set dep [string_change $dep $changelist]
		}
		set list [regexp -all -inline {\$\{?@[a-zA-z0-9_]+\(?\)?\}?} $dep]
		foreach pattern $list {
			regexp {[a-zA-z0-9_]+} $pattern base
			if {![info exists ::$base]} {
				cgmaketarget @$base newids
				cgmakewait $newids
				if {![info exists ::@$base]} {
					exiterror "Could not create variable @$base"
				}
			}
			set values [get ::@$base]
			set newdep {}
			foreach e $dep {
				foreach v $values {
					lappend newdep [string_change $e [list $pattern $v]]
				}
			}
			set dep $newdep
		}
		lappend result {*}$dep
#		foreach e $dep {
#			if {[string first * $e] != -1} {
#				set g [glob -nocomplain $e]
#				if {[llength $g]} {
#					lappend result {*}$g
#				} else {
#					lappend result $e
#				}
#			} else {
#				lappend result $e
#			}
#		}
	}
	return $result
}

proc cgmakewait {ids} {
}

proc cgmaketarget {target newidsVar} {
	global cgmakedata cgmakeargs cgmakeids cgmakeroot
	upvar $newidsVar newids
	cd $cgmakeroot
	set newids {}
	set target [uplevel #0 [list subst -nobackslashes -nocommands $target]]
	if {[llength [glob -nocomplain $target]]} {
		putslog "file $target exists"
		return ok
	} elseif {[info exists cgmakeids($target)]} {
		if {$cgmakeids($target) eq ""} {
			putslog "Target $target already done"
		} else {
			putslog "Target $target already submitted (id $cgmakeids($target))"
			lappend newids $cgmakeids($target)
		}
		return ok
	}
	putslog "---------- Finding rules for target $target ----------"
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
		putslog "cannot make $target"
		return "cannot make $target"
	}
	list_foreach {pattern deps vars precode code submitopts} $checklist {
		putslog "----- Test rule $pattern <-- $deps -----"
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
		set deps [uplevel #0 [list subst -nobackslashes -nocommands $deps]]
		set ids {}
		set depsok 1
		set deps [cgmake_expanddeps $deps]
		if {$precode ne ""} {
			putslog "Running precode: $precode"
			set precode [cgmakecode $precode $target $targetvars $deps $vars]
			proc cgmake_temp {} $precode
			uplevel #0 cgmake_temp
		}
		foreach dep $deps {
			set status [cgmaketarget $dep newids]
			if {[regexp {^cannot make} $status]} {
				set depsok 0
				break
			}
			if {[llength $newids]} {
				lappend ids {*}$newids
			}
		}
		if {$depsok} break
	}
	putslog "-------------------- Target $target --------------------"
	set ecode [cgmakecode $code $target $targetvars $deps $vars]
	if {$cgmakeargs(distribute) == 0} {
		incr cgmakeids()
		set cgmakeids($target) $cgmakeids()
		proc cgmake_temp {} $ecode
		putslog "Making $target"
		uplevel #0 cgmake_temp
		set cgmakeids($target) {}
		putslog "Made $target\n"
	} else {
		exiterror "other than -d 0 not implemented yet"
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
	global cgmakedata cgmakeargs
	unset -nocomplain cgmakedata
	unset -nocomplain cgmakeargs
}

proc cgmake {args} {
	global cgmakeargs cgmakeids cgmakeroot
	unset -nocomplain cgmakeids
	set cgmakeroot [pwd]
	set cgmakeids() 0
	set args [cgmake_args {*}$args]
	if {![llength $args]} {set args all}
	foreach target $args {
		cgmaketarget $target newids
	}
}

if 0 {
	set cgmakefile ../cgmake.tcl ../tests/data test
}

