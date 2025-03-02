# create soft link (dest points to src) using relative path (unless absolute == 1)
# allow links to non-existing files
proc mklink {args} {
	set absolute 0
	set matchtime 1
	cg_options mklink args {
		-absolute {
			set absolute $value
		}
		-matchtime {
			set matchtime $value
		}
	} {src dest absolute} 2 3 {
		make a soflink (dest points to src)
	}
	set src [file_absolute $src]
	set keepsrc $src
	set dest [file_absolute $dest]
	if {$src eq $dest} {
		error "cannot mklink a file to itself: $src"
	}
	if {!$absolute} {
		set pos 0
		set ssrc [file split $keepsrc]
		set sdest [file split $dest]
		# puts $ssrc\n$sdest
		foreach s $ssrc d $sdest {
			if {$s ne $d} break
			incr pos
		}
		if {$s eq "" && $pos != 0} {incr pos -1}
		if {$pos > 1} {
			set prelen [expr {[llength $sdest]-$pos -1}]
			set src [file join {*}[list_fill $prelen ..] {*}[lrange $ssrc $pos end]]
		}
	}
	# cannot do file exists first because a (existing) link to an unexisting file would return 0
	set err [catch {file link $dest} link]
	if {$err} {
		if {[file exists $dest]} {
			error "cg mklink error: destination exists and is not a link: $dest"
		}
		set make 1
	} elseif {$link ne "$src"} {
		set make 1
	} else {
		set make 0
	}
	if {$make} {
		file delete $dest
		if {[file exists $keepsrc]} {
			file link -symbolic $dest $src
		} else {
			set keeppwd [pwd]
			mkdir [file dir $dest]
			cd [file dir $dest]
			# for some reason, linking to a non-existing file this way (called from Tcl)
			# sometimes gives an "Error: couldn't execute "ln": too many levels of symbolic links" error
			# workaround: called from bash -> no error
			if {[catch {
				exec ln -sf $src [file tail $dest]
			}]} {
				exec bash -c "ln -sf \'$src\' \'[file tail $dest]\'"
			}
			cd $keeppwd
		}
	}
	if {[file exists $keepsrc] && $matchtime && ![catch {file lstat $keepsrc a}]} {
		set mtime $a(mtime)
		if {[file_mtime $dest] ne $mtime} {
			exec touch -h -d [clock format $mtime] $dest
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

proc mklink_asjob {dep target} {
	upvar job_logdir job_logdir
	job mklink-$target -checkcompressed 0 -deps {$dep} -targets {$target} -code {
		mklink $dep $target
	}
}

proc cg_mklink {args} {
	set absolute 0
	set matchtime 1
	cg_options mklink args {
		-absolute - -a {
			set absolute $value
		}
		-matchtime {
			set matchtime $value
		}
	} {src dest} 2 ... {
		make a soflink (dest points to src)
	}
	if {[llength $args]} {
		set srcs [list $src $dest]
		set dest [lindex $args end]
		lappend srcs {*}[lrange $args 0 end-1]
		if {![file isdir $dest]} {
			error "only can make links from multiple files in a directory"
		}
		foreach src $srcs {
			mklink -absolute $absolute -matchtime $matchtime $src $dest/[file tail $src]
		}
	} else {
		mklink -absolute $absolute -matchtime $matchtime $src $dest
	}
}
