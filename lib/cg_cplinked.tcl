# rmlinked removes a directory, but only if it contains only (soft)links, or directories containing softlinks
proc rmlinked {src} {
	set src [file join [pwd] $src]
	if {![file isdir $src]} {
		if {![catch {file link $src} link]} {
			file delete $src
		}
		return
	}
	# exec chmod g+w $dest
	set files [glob -nocomplain $src/*]
	set keep 0
	foreach file $files {
		if {[file isdir $file] && [catch {file link $file} link]} {
			if {[rmlinked $file]} {set keep 1}
		} else {
			if {![catch {file link $file} link]} {
				file delete $file
			} else {
				set keep 1
			}
		}
	}
	if {!$keep} {file delete $src}
	return $keep
}

proc cg_rmlinked {args} {
	if {[llength $args] < 1} {
		exiterror {wrong # args: should be "cg rmlinked dir ..."}
	}
	foreach {src} $args {
		rmlinked $src
	}
}

# cplinked copies a directory, but by making (soft)links to the original directory for all files
# existing files (that are not links) will be backed up
proc cplinked {src dest {absolute 0}} {
	set src [file join [pwd] $src]
	set dest [file join [pwd] $dest]
	if {![file isdir $src] || ![catch {file link $src} link]} {
		cplinked_file $src $dest $absolute
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
			cplinked $file $destfile $absolute
		} else {
			cplinked_file $file $destfile $absolute
		}
	}
}

proc cplinked_file {file destfile {absolute 0}} {
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
			file rename -force $destfile $renamefile$num
		}
	}
	if {!$absolute} {
		mklink $file $destfile
	} else {
		file link -symbolic $destfile [file_absolute $file]
	}
}

proc cg_cplinked {args} {
	if {[lindex $args 0] eq "-absolute"} {
		set absolute 1
		set args [lrange $args 1 end]
	} else {
		set absolute 0
	}
	if {[llength $args] != 2} {
		exiterror {wrong # args: should be "cg cplinked src dest"}
	}
	foreach {src dest} $args break
	if {[file isdir $dest]} {
		set tail [file tail $src]
		set dest $dest/$tail
	}
	cplinked $src $dest $absolute
}
