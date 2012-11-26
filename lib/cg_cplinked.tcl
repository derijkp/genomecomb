# cplinked copies a directory, but by making (soft)links to the original directory for all files
# existing files (that are not links) will be backed up
proc cplinked {src dest} {
	set src [file normalize $src]
	set dest [file normalize $dest]
	if {![file isdir $src]} {
		cplinked_file $src $dest
		return
	}
	if {[file exists $dest]} {
		if {![file isdir $dest]} {error "cannot overwrite non-directory `$dest' with directory '$src'"}
	}
	file mkdir $dest
	# exec chmod g+w $dest
	set files [glob -nocomplain $src/*]
	foreach file $files {
		set destfile $dest/[file tail $file]
		if {[file isdir $file]} {
			cplinked $file $destfile
		} else {
			cplinked_file $file $destfile
		}
	}
}

proc cplinked_file {file destfile} {
	puts "$file -> $destfile"
	if {[file exists $destfile]} {
		if {![catch {file link $destfile} link]} {
			file delete $destfile
		} else {
			set renamefile $destfile.old
			set num 1
			while {[file exists $renamefile$num]} {
				incr num
			}
			puts "destfile $destfile exists, renamed to $renamefile$num"
			file rename $destfile $renamefile$num
		}
	}
	mklink $file $destfile
}

proc cg_cplinked {args} {
	if {[llength $args] != 2} {
		exiterror {wrong # args: should be "cg cplinked src dest"}
	}
	foreach {src dest} $args break
	if {[file isdir $dest]} {
		set tail [file tail $src]
		set dest $dest/$tail
	}
	cplinked $src $dest
}
