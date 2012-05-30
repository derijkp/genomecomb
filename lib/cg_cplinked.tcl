# cplinked copies a directory, but by making (soft)links to the original directory for all files
# existing files (that are not links) will be backed up
proc cplinked {src dest} {
	if {[file pathtype $dest] ne "absolute"} {
		set dest [file normalize $dest]
	}
	file mkdir $dest
	# exec chmod g+w $dest
	set keeppwd [pwd]
	if {[file pathtype $src] ne "absolute"} {
		set src ../$src
		cd $dest
	}
#	# find common dir (later)
#	set num 0
#	set commondir {}
#	foreach s [lrange [file split $src] 1 end] d [lrange [file split $dest] 1 end] {
#		if {$s ne $d} break
#		set commondir $commondir/$s
#		incr num
#	}
	set files [glob -nocomplain $src/*]
	foreach file $files {
		set destfile $dest/[file tail $file]
		if {[file isdir $file]} {
			cplinked $file $destfile
		} else {
			if {[file exists $destfile]} {
				if {![catch {file link $destfile} link]} {
					file delete $destfile
				} else {
					set renamefile $destfile.old
					set num 1
					while {[file exists $renamefile$num]} {
						incr num
					}
					puts stderr "destfile $destfile exists, renamed to $renamefile$num"
					file rename $destfile $renamefile$num
				}
			}
			file link $destfile $file
			# exec ln -s $file $destfile
		}
	}
	cd $keeppwd
}

proc cg_cplinked {args} {
	if {[llength $args] != 2} {
		exiterror {wrong # args: should be "cg cplinked src dest"}
	}
	foreach {src dest} $args break
	cplinked $src $dest
}
