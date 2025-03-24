proc file_absolute {file} {
	if {$file eq ""} {return $file}
	if {[string index $file 0] eq "-" && [file_root $file] eq "-"} {
		return $file
	}
	if {[string index $file 0] eq "~"} {
		if {[string length $file] == 1} {
			return $::env(HOME)
		} elseif {[string index $file 1] eq "/"} {
			set file [file join $::env(HOME) [string range $file 2 end]]
		} else {
			set file ./$file
		}
	}
	set result {}
	foreach el [file split [file join [pwd] $file]] {
		if {$el eq ".."} {
			if {[llength $result] <= 1} {error "file_absolute error: cannot .. past root"}
			list_pop result
		} elseif {$el ne "." && $el ne ""} {
			lappend result $el
		}
	}
	file join {*}$result
}

proc file_sample {file} {
	join [lrange [split [file root [file tail [gzroot $file]]] -] end end] -
}

proc file_analysis {file} {
	join [lrange [split [file root [file tail [gzroot $file]]] -] 1 end] -
}

namespace eval genomecomb {}
if {[info commands genomecomb::cd.ori] eq ""} {
	rename cd genomecomb::cd.ori
}
if {[info commands genomecomb::pwd.ori] eq ""} {
	rename pwd genomecomb::pwd.ori
}
proc cd {path} {
	set path [file_absolute $path]
	genomecomb::cd.ori $path
	set ::genomecomb::cwd $path
	set ::env(PWD) $path
}

proc pwd {} {
	if {![info exists ::genomecomb::cwd]} {
		if {[info exists ::env(PWD)]} {
			set ::genomecomb::cwd $::env(PWD)
		} else {
			set ::genomecomb::cwd [::genomecomb::pwd.ori]
		}
	}
	return $::genomecomb::cwd
}

proc mkdir {dir} {
	if {![file exists $dir]} {file mkdir $dir}	
}

proc rm {args} {
	set opts {}
	set pos 0
	set force 0
	set warning 0
	set recursive 0
	cg_options rm args {
		-force {set force $value}
		-recursive {set recursive $value}
		-warning {set warning $value}
	} file 1 ... {
		delete files, includes support for deleting shadow dirs
	}
	set files [list $file {*}$args]
	foreach file $files {
		if {[file isdir $file] && ![catch {file link $file} shadow] && [file exists $shadow/shadow_source]} {
			file delete -force $shadow
			if {[catch {file delete -force $shadow} msg]} {
				if {![regexp "no such file or directory" $msg]} {
					if {$warning} {
						putswarning "warning: could not delete shadow dir $shadow of file $file: $msg"
					} else {
						error "error deleting shadow dir $shadow of file $file: $msg"
					}
				}
			}
		}
		if {$recursive && [file isdir $file]} {
			foreach subfile [glob -nocomplain $file/*] {
				rm -force $force -recursive $recursive -warning $warning $subfile
			}
		}
		if {[catch {file delete {*}$opts $file} msg]} {
			if {![regexp "no such file or directory" $msg]} {
				if {$warning} {
					putswarning "warning: could not delete file $file: $msg"
				} else {
					error "error deleting file $file: $msg"
`				}
			}
		}
	}
}

proc file_mtime {file} {
	file lstat $file a
	return $a(mtime)
}

proc pathsep {} {
	if {$::tcl_platform(platform) eq "windows"} {return \;} else {return \:}
}

proc copywithindex {file target args} {
	file copy -force $file $target
	foreach {file2 target2} $args {
		if {[file exists $file2]} {
			file copy -force $file2 $target2
		} elseif {[file extension $target2] eq ".tbi" && [file extension [file root $target2]] eq ".gz"} {
			exec tabix -p vcf [file root $target2]
		}
	}	
}

proc file_part {file {start {}} {end {}}} {
	if {$start eq ""} {
		set start end ; set end end
	} elseif {$end eq ""} {
		set end $start
	}
	file join {*}[lrange [file split $file] $start $end]
}

proc workdir {file} {
	set workdir [gzroot $file].temp
	if {[file exists $workdir] && ![file isdir $workdir]} {
		file delete $workdir
	}
	file mkdir $workdir
	job_cleanup_add $workdir
	return $workdir
}

proc shadow_workdir {file} {
	set workdir [gzroot $file].temp
	if {[file exists $workdir] && ![file isdir $workdir]} {
		file delete $workdir
	}
	shadow_mkdir $workdir
	job_cleanup_add_shadow $workdir
	return $workdir
}

proc touch {args} {
	exec touch {*}$args
}

# like file copy, but always used -force, and if target is a file, it changes it's time to the current
proc file_copy {args} {
	file copy -force {*}$args
	set target [lindex $args end]
	exec touch $target
}

proc follow_links {file} {
	set file [file_absolute $file]
	while 1 {
		if {[catch {set file [file join [file dir $file] [file readlink $file]]}]} break
	}
	return $file
}

proc sclose {f} {
	if {$f in "stdin stdout stderr"} return
	close $f
}

proc read_tsv {file {fields {}}} {
	set f [gzopen $file]
	set header [tsv_open $f]
	if {$fields eq ""} {
		set result [list $header]
		while {[gets $f line] != -1} {
			lappend result [split $line \t]
		}
	} else {
		set poss [list_cor $header $fields]
		set result {}
		while {[gets $f line] != -1} {
			lappend result [list_sub [split $line \t] $poss]
		}
	}
	gzclose $f
	return $result
}

proc file_or_link_exists {file} {
	if {[file exists $file]} {return 1}
	if {[catch {file link $file}]} {return 0} else {return 1}
}

proc file_write_gz {args} {
	set list [list_pop args]
	set filename [list_pop args]
	set f [wgzopen $filename]
	eval fconfigure $f -buffersize 100000 $args
	puts -nonewline $f $list
	close $f
}

if {[info commands genomecomb::open.ori] eq ""} {
	rename open genomecomb::open.ori
}
proc open args {
	set f [genomecomb::open.ori {*}$args]
	set ::genomecomb_gzopen_info($f) $args
	return $f
}

if {[info commands genomecomb::close.ori] eq ""} {
	rename close genomecomb::close.ori
}
proc close f {
	if {$f in "stdin stdout stderr"} return
	genomecomb::close.ori $f
	unset -nocomplain ::genomecomb_gzopen_info($f)
}
