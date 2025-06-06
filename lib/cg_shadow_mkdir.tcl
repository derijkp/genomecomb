proc shadow_clean {{shadowdir {}}} {
	global env
	if {$shadowdir ne ""} {
		# shadowdir given directly, do nothing extra
	} elseif {[info exists env(SHADOWDIR)]} {
		set shadowdir $env(SHADOWDIR)
	} else {
		error "shadowdir not given, and could not find env var SHADOWDIR"
	}
	set shadowdir [file_absolute $shadowdir]
	foreach shadow [glob $shadowdir/*] {
		# if link exists, but points to non-existing file, file exists wil return 0
		# so check with file link as well for that
		if {[catch {file link $shadow/shadow_source} link] && ![file exists $shadow/shadow_source]} {
			puts "\nskipping $shadow : does not have shadow_source"
			continue
		}
		if {
			[file exists $link] &&
			![catch {file link $link} linklink] &&
			$linklink eq $shadow
		} continue
		puts -nonewline .
		file delete -force $shadow
	}
}

proc cg_shadow_clean {args} {
	set shadowdir {}
	cg_options shadow_clean args {
		-shadowdir {
			set shadowdir $value
		}
	} {shadowdir} 0 1 {
		clean up shadowdirs that are no longer used
	}
	shadow_clean $shadowdir 
}

# shadow_delete is folded into rm, so no longer really needed
# if warning = 1, a warning is sent to stderr if the file cannot be deleted, but no error is generated
proc shadow_delete {link {warning 0}} {
	if {![file isdir $link]} {
		rm -warning $warning $link
		return
	}
	if {[catch {file link $link} shadow]} {
		rm -warning $warning -force 1 -recursive 1 $link
		return
	}
	if {![file exists $shadow/shadow_source]} {
		puts stderr "$link links to $shadow is not a shadow (file $shadow/shadow_source does not exist), only removing link, not target"
	} else {
		rm -warning $warning -force 1 -recursive 1 $shadow
	}
	rm -warning $warning $link
}

proc cg_shadow_delete {args} {
	cg_options shadow_delete args {
	} {link} 1 ... {
		delete shadow dirs
	}
	rm -force 1 -recursive 1 $link {*}$args
}

# shadow_delete is folded into rm ->
# cg_rm is mostly the same as cg_shadow_delete (cg_shadow_delete uses -force by default)
proc cg_rm {args} {
	rm {*}$args
}

# shadow_mkdir creates a temp or workdir for distributed analysis.
# - it can be accessed from different nodes in cluster
# - it has to be explicitely deleted: it does not get auto-deleted when a job or program finishes
proc shadow_mkdir {link {shadowdir {}}} {
	global env
	if {[file exists $link]} {
		return $link
	}
	if {![file exists [file dir $link]]} {
		file mkdir [file dir $link]
	}
	if {$shadowdir ne ""} {
		# shadowdir given directly, do nothing extra
	} elseif {[info exists env(SHADOWDIR)]} {
		set shadowdir $env(SHADOWDIR)
	} else {
		# putslog "shadowdir not given, and could not find evn var SHADOWDIR, using $link directly"
		mkdir $link
		return $link
	}
	set shadowdir [file_absolute $shadowdir]
	# find unused dir for shadow and create
	for {set i 0} {$i < 20} {incr i} {
		set shadow [file join $shadowdir shadow.[pid]-[Extral::randstring 20]]-[file tail $link]
		if {[file exists $shadow]} continue
		if {[catch {
			file mkdir $shadow
			if {$::tcl_platform(platform) eq "unix"} {
				file attributes $shadow -permissions 0700
			}
			set files [glob -nocomplain $shadow/*]
			if {[llength $files]} {
				error "Very fishy: there are files in the temporary directory I just created"
			}
		}]} continue
		break
	}
	# create backlink (for cleaning)
	mklink $link $shadow/shadow_source 1
	# create link
	mklink $shadow $link 1
	return $link
}

proc cg_shadow_mkdir {args} {
	set shadowdir {}
	cg_options shadow_mkdir args {
		-shadowdir {
			set shadowdir $value
		}
	} {link shadowdir} 1 2 {
		make a shadow dir
	}
	shadow_mkdir $link $shadowdir
}
