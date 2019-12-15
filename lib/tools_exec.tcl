proc chanexec {in out pipe} {
	if {$pipe eq ""} {
		set o $out
	} else {
		set o [open [list | {*}$pipe >@ $out 2>@ stderr] w]
	}
	if {[info exists ::filebuffer($in)]} {
		foreach line $::filebuffer($in) {
			puts $o $line
		}
		unset ::filebuffer($in)
	}
	catch_fileaccess {fcopy $in $o} $in $o
	if {$in ne "stdin"} {gzclose $in}
	# gzclose would catch error "child process exited abnormally", which would cause us to miss errors in the pipe
	# we only want to do this for input (interupted decompression)
	# gzcloseout only removes the "child process exited abnormally" from the error, but propagates the exit error
	if {$o ne "stdout"} {gzcloseout $o}
}

proc catch_exec {args} {
	set error [catch {
		exec {*}$args
	} msg opt]
	if {$error} {
		if {$::errorCode ne "NONE"} {
			dict unset opt -level
			set errorInfo "$msg\n    while executing\n$args"
			return -code $error -errorcode $::errorCode -errorinfo $errorInfo $msg
		}
	}
	return $msg
}

proc catch_fileaccess {cmd args} {
	set error [catch {uplevel $cmd} msg]
	if {$error} {
		set keeperrorCode $::errorCode
		set keeperrorInfo $::errorInfo
		foreach f $args {
			if {$f in "stdin stdout stderr"} continue
			if {[catch {gzclose $f} msg]} {
				append error "\nerror closing $f: $msg"
			}
		}
		set ::errorCode $keeperrorCode
		set ::errorInfo $keeperrorInfo
		error $msg
	}
}

proc catchchildkilled_exec {args} {
	if {[catch {
		exec {*}$args
	} msg opt]} {
		if {$::errorCode ne "NONE" && ![string match {CHILDKILLED * SIGPIPE *} $::errorCode]} {
			dict unset opt -level
			return -options $opt $msg
		}
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
		set bindir $appdir/bin
	} elseif {[file tail $appdir] eq "cg_viz"} {
		# we are being run from dev cg_viz
		set genomecombdir [file dir $appdir]
		set externdir $genomecombdir/extern
		set bindir $genomecombdir/bin
	} else {
		# we are being run from dev 
		set genomecombdir $appdir
		set externdir $genomecombdir/extern
		set bindir $genomecombdir/bin
	}
	set env(PATH) $bindir[pathsep]$externdir[pathsep]$genomecombdir[pathsep]$env(PATH)
	if {[info exists env(LD_LIBRARY_PATH)]} {
		set env(LD_LIBRARY_PATH) $::externdir/lib:$env(LD_LIBRARY_PATH)
	} else {
		set env(LD_LIBRARY_PATH) $::externdir/lib
	}
	return $genomecombdir
}

