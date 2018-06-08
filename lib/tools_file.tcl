#doc scratchdir title {
#scratchdir
#} shortdescr {
# returns a directory in which temporary files can be stored. This directory is specific to one proces:
# no other processes will (should) write in this directory. Subsequent calls to the function within one process
# will allways be the same directory, The program has to take care not to overwrite its own files
# Temporary files returned by tempfile are also in this directory (named like _Extral_temp_1.tmp)
# The program should also not overwrite these
# The temporary directory is deleted when the program exits by an atexit handler
# 
#}

proc scratchdir {} {
	global env
	if {![info exists env(SCRATCHDIR)]} {
		# putslog "Could not find SCRATCHDIR, using tempdir"
		set tempdir [tempdir]
		if {[file isdir [file dir $tempdir]/scratch]} {
			set env(SCRATCHDIR) [file dir $tempdir]/scratch
		} else {
			return $tempdir
		}
	}
	if {![info exists ::Extral::scratchdir]} {
		for {set i 0} {$i < 20} {incr i} {
			set testdir [file join $env(SCRATCHDIR) scratchExtral.[pid]-[Extral::randstring 20]]
			if {[file exists $testdir]} continue
			if {[catch {
				file mkdir $testdir
				if {$::tcl_platform(platform) eq "unix"} {
					file attributes $testdir -permissions 0700
				}
				set files [glob -nocomplain $testdir/*]
				if {[llength $files]} {
					error "Very fishy: there are files in the temporary directory I just created"
				}
				set ::Extral::scratchdir $testdir
				set ::Extral::scratchnum 0
			}]} continue
			break
		}
	}
	if {![info exists ::Extral::scratchdir]} {
		error "couldn't create scratch directory in $env(SCRATCHDIR) (defined by env variable SCRATCHDIR)"
	}
	return $::Extral::scratchdir
}

proc scratchfile {{action {get}} {type file}} {
	switch $action {
		get {
			set scratchdir [scratchdir]
			return [file join $scratchdir _Extral_scratch_[incr ::Extral::scratchnum].tmp]
		}
		clean {
			set scratchdir [scratchdir]
			catch {file delete -force $scratchdir}
			unset ::Extral::scratchdir
		}
		cleanall {
			catch {file delete -force {*}[glob [file join $scratch_dir scratchExtral*]]}
		}
		default {
			return -code error "bad option \"$action\": must be get, clean or cleanall"
		}
	}
}

proc filetemp {file {write 1}} {
	if {![file exists $file.temp]} {
		set result $file.temp
	} else {
		set num 2
		while {[file exists $file.temp$num]} {incr num}
		set result $file.temp$num
	}
	if {$write} {file_write $result {}}
	return $result
}

proc filetemp_ext {file {write 1}} {
	set ext [file extension $file]
	if {![isgzext $ext]} {set ext {}}
	if {![file exists $file.temp$ext]} {
		set result $file.temp$ext
	} else {
		set num 2
		while {[file exists $file.temp$num$ext]} {incr num}
		set result $file.temp$num$ext
	}
	if {$write} {file_write $result {}}
	return $result
}

proc maxopenfiles {{force 0}} {
	global maxopenfiles
	if {$force} {unset -nocomplain maxopenfiles}
	if {[info exists maxopenfiles]} {
		if {$force} {
			unset -nocomplain maxopenfiles
		} else {
			return $maxopenfiles
		}
	}
	if {[file exists /proc/self/limits]} {
		set c [file_read /proc/self/limits]
		if {[regexp {Max open files  *([0-9]+)} $c temp maxopenfiles]} {
			incr maxopenfiles -10
			return [max $maxopenfiles 10]
		}
	}
	if {![catch {exec sh -c {ulimit -n}} temp] && [isint $temp]} {
		set maxopenfiles [expr {$temp - 10}]
		return [max $maxopenfiles 10]
	}
	return 1000
}

proc file_add {file args} {
	set o [open $file a]
	foreach arg $args {
		puts $o $arg
	}
	close $o
}

proc checkfile {args} {
	foreach filename $args {
		if {![catch {glob $filename} list]} {
			return [lindex $list 0]
		}
	}
	return [lindex $args 0]
}

proc checkfiles {args} {
	set result {}
	foreach filename $args {
		if {![catch {glob $filename} list]} {
			lappend result {*}$list
		}
	}
	return $result
}

proc getline f {
	set line [split [gets $f] \t]
	while {![llength $line] && ![eof $f]} {
		set line [split [gets $f] \t]
	}
	return $line
}

proc getnotempty {f} {
	set line {}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line]} break
	}
	return $line
}

# returns absolute path of waht file is linking to
proc file_resolve {file {lastlinkVar {}}} {
	if {$lastlinkVar ne ""} {upvar $lastlinkVar lastlink}
	if {$::tcl_platform(platform) eq "unix"} {
		set file [file_absolute $file]
		while 1 {
			if {[catch {set link [file readlink $file]}]} break
			if {[file pathtype $link] ne "absolute"} {set link [file_absolute [file join [file dir $file] $link]]}
			set lastlink $file
			set file [file_absolute $link]
		}
	}
	return $file
}

# returns absolute path of what file is linking to
proc find_link {file {level {}}} {
	set file [file_absolute $file]
	while 1 {
		if {[catch {
			set file [file join [file dir $file] [file readlink $file]]
			set file [file_absolute $file]
		}]} break
		if {$level ne "" && ![incr level -1]} break
	}
	return $file
}

# make softlink as in tcl (linkname (1st) points to linkdest),
# but linkdest is adjusted relative to linkname, and does not have to exist
proc file_link {linkname linkdest} {
	set dir [file dir $linkname]
	if {[file pathtype $linkdest] eq "absolute"} {
		set dest $linkdest
	} else {
		set dest [file join $dir $linkdest]
	}
	if {![file exists $dest]} {
		file_write $dest temp
		file link $linkname $linkdest
		file delete $dest
	} else {
		file link $linkname $linkdest
	}
}

# copy recursively with permissions and dates using hardlinks
proc hardlink {args} {
	if {[llength $args] < 2} {error "wrong # args: should be \"hardlink src ... dest\""}
	exec cp -al {*}$args
}

# copy recursively with permissions and dates using hardlinks, but if that does not work, copy normally
proc hardcopy {args} {
	if {[llength $args] < 2} {error "wrong # args: should be \"hardcopy src ... dest\""}
	if {[catch {hardlink {*}$args}]} {
		exec cp -a {*}$args
	}
}

# create soft link (dest points to src) using relative path (unless absolute == 1)
# allow links to non-exusting files
proc mklink {src dest {absolute 0}} {
	set src [file_absolute $src]
	set dest [file_absolute $dest]
	if {!$absolute} {
		set pos 0
		set ssrc [file split $src]
		set sdest [file split $dest]
		# puts $ssrc\n$sdest
		foreach s $ssrc d $sdest {
			if {$s ne $d} break
			incr pos
		}
		if {$pos > 1} {
			set prelen [expr {[llength $sdest]-$pos -1}]
			set src [file join {*}[list_fill $prelen ..] {*}[lrange $ssrc $pos end]]
		}
	}
	set err [catch {file link $dest} link]
	if {!$err || $link ne "$src"} {
		file delete $dest
	}
	if {![file exists $dest]} {
		if {[file exists $src]} {
			file link -symbolic $dest $src
		} else {
			set keeppwd [pwd]
			cd [file dir $dest]
			exec ln -s $src [file tail $dest]
			cd $keeppwd
		}
	}
}

# same as mklink, but add proper extension to dest if src is compressed (and dest does not have already the correct compression extension)
proc gzmklink {src dest} {
	set src [gzfile $src]
	set ext_s [file extension $src]
	set ext_d [file extension $dest]
	if {$ext_s ne $ext_d && [isgzext $ext_s]} {
		mklink $src $dest$ext_s
	} else {
		mklink $src $dest
	}
}

proc mklink_job {dep target} {
	upvar job_logdir job_logdir
	job mklink-$target -checkcompressed 0 -deps {$dep} -targets {$target} -code {
		mklink $dep $target
	}
}

proc file_timestamp {file} {
	clock format [file mtime $file] -format "%Y-%m-%d_%H_%M_%S"
}

proc getlink {file} {
	file_absolute [file join [file dir $file] [file link $file]]
}

proc gzlink {file} {
	if {[file exists $file]} {return $file}
	if {[catch {glob $file}]} {return ""}
	gzfile [getlink $file]
}

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

proc filename {file} {
	if {![catch {file link $file} link]} {
		return [file tail $link]
	} else {
		return [file tail $file]
	}
}

proc file2refname {file} {
	if {![catch {file link $file} link]} {
		set tail [file tail $link]
	} else {
		set tail [file tail $file]
	}
	regexp {genome_([^_.]+)\.ifas} $tail temp tail
	return $tail
}

proc file_rootname {file} {
	set base [file root [gzroot [file tail $file]]]
	set root [join [lrange [split $base -] 1 end] -]
	if {$root ne ""} {return $root}
	return $base
}
