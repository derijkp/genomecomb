proc chanexec {in out pipe} {
	set o [open "|\ $pipe\ >@\ $out 2>@\ stderr" w]
	if {[info exists ::filebuffer($in)]} {
		foreach line $::filebuffer($in) {
			puts $o $line
		}
		unset ::filebuffer($in)
	}
	fcopy $in $o
	if {$in ne "stdin"} {catch {gzclose $in}}
	if {$out ne "stdout"} {catch {close $out}}
	close $o
}

proc catchstderr_exec {args} {
	if {[catch {
		exec {*}$args 2>@1
	} msg]} {
		error $msg
	}
	return $msg
}

proc progress {cmd args} {
	if {![llength [info commands winfo]]} {
		global progresslevel
		if {![info exists progresslevel]} {
			set progresslevel -1
		}
		# no Tk, so we are running commandline: ignore a lot of the progress commands
		switch $cmd {
			start {
				incr progresslevel
			}
			stop {
				if {$progresslevel == -1} return
				incr progresslevel -1
			}
			protect {
				foreach {code on_error} $args break
				set error [catch {uplevel 1 $code} result]
				if {$error} {
					if {$on_error eq ""} {
						set errorInfo $::errorInfo
						return -code error -errorinfo $errorInfo $result
					} else {
						set ::errorResult $result
						uplevel 1 $on_error
					}
				} else {
					return $result
				}
			}
		}
	} else {
		switch $cmd {
			onerror {
				Classy::Progress on_error [subst {
					[lindex $args 0]
				}]
			}
			startdisplay {
				if {([Classy::Progress level] == -1)} {
				}
				uplevel 1 [list Classy::Progress start] $args
			}
			start {
				if {([Classy::Progress level] == -1)} {
				}
				uplevel 1 [list Classy::Progress $cmd] $args
			}
			stop {
				Classy::Progress stop
				if {([Classy::Progress level] == -1)} {
					Classy::Progress on_error {}
				}
			}
			next {
				uplevel 1 [list Classy::Progress $cmd] $args
				update idletasks
			}
			cancel {
				Classy::Progress cancel
			}
			default {
				uplevel 1 [list Classy::Progress $cmd] $args
			}
		}
	}
}

proc exiterror errormessage {
	puts stderr $errormessage
	exit 1
}

proc reload {} {
	global appdir
	foreach file [glob $appdir/lib/*.tcl $appdir/lib-exp/*.tcl] {
		puts "sourcing $file"
		uplevel 0 source $file
	}
}

proc genomecombenv {} {
	global auto_path env appdir tcl_dirtcl genomecombdir externdir
	if {![info exists appdir]} {
		set appdir ${genomecomb::dir}
	}
	if {[file dir [file dir $appdir]] eq [get tcl_dirtcl ""]} {
		# we are being run from a dirtcl installation in apps/cg
		set genomecombdir $tcl_dirtcl
		set externdir $genomecombdir/bin
		set env(PATH) $appdir/bin:$externdir:$genomecombdir:$env(PATH)
	} elseif {[file tail $appdir] eq "cg_viz"} {
		# we are being run from dev cg_viz
		set genomecombdir [file dir $appdir]
		set externdir $genomecombdir/extern
		set env(PATH) $genomecombdir/bin:$externdir:$genomecombdir:$env(PATH)
	} else {
		# we are being run from dev 
		set genomecombdir $appdir
		set externdir $genomecombdir/extern
		set env(PATH) $genomecombdir/bin:$externdir:$genomecombdir:$env(PATH)
	}
	return $genomecombdir
}
