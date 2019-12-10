proc infofile_read {file} {
	set f [open $file]
	set header [tsv_open $f]
	if {$header ne {key value}} {
		error "file $file has a wrong header for an info file, should be: key value"
	}
	set result {}
	while {![eof $f]} {
		set line [gets $f]
		if {[regexp {^([^\t]+)\t(.*)$} $line temp key value]} {
			dict set result $key $value
		}
	}
	return $result
}

proc infofile_write {file data} {
	set o [open $file.temp w]
	puts $o key\tvalue
	foreach {key value} $data {
		puts $o $key\t$value
	}
	close $o
	file rename -force -- $file.temp $file
}

# initialize variables (from) projectinfo.tsv file
# dir is the project directory. The command will search for the file projectinfo.tsv in it. If it does not exist, it is created.
# Further arguments a variables, that will be assigned values after projectinfo:
# If a variable already has a value before projectinfo
#   If the variable is already defined in projectinfo.tsv, and differs from the given value, an error is given
#   If the variable is not in projectinfo.tsv yet, it is stored there
# If the variable is not defined before projectinfo, it should have a default value (list of variable name and default value)
#   It it is in projectinfo.tsv, the value in projectinfo.tsv will be assigned
#   It it is not in projectinfo.tsv, the default value will be assigned and stored in projectinfo.tsv
proc projectinfo {dir args} {
	set projectinfofile $dir/projectinfo.tsv
	# check projectinfo
	if {[file exists $projectinfofile]} {
		set infod [infofile_read $projectinfofile]
	} else {
		set infod {}
	}
	foreach line $args {
		set def {}
		foreach {varVar def} $line break
		if {[llength $line] > 1} {set usedef 1} else {set usedef 0}
		upvar $varVar var
		if {![info exists var] || $var eq {}} {
			if {![dict exists $infod $varVar]} {
				if {$usedef} {
					set var $def
				} else {
					error "error: no $varVar parameter given, and it is also not in the projectinfo.tsv"
				}
			} else {
				set var [dict get $infod $varVar]
			}
		} else {
			if {[dict exists $infod $varVar] && $var ne [dict get $infod $varVar]} {
				error "error: The $varVar parameter given ($var) differs from the one in the projectinfo file $projectinfofile ([dict get $infod $varVar])"
			}
		}
		dict set infod $varVar $var
	}
	infofile_write $projectinfofile $infod
	return $infod
}

proc testmultitarget {target files} {
	file delete $target.temp
	if {[file exists $target]} {
		# test if existing target is already ok
		set ok 1
		set done [cg select -a $target]
		foreach testfile $files {
			if {![file exists $testfile] || [file mtime $target] < [file mtime $testfile]} {
				putslog "$testfile is newer than $target"
				set ok 0
			}
			set name [file_rootname $testfile]
			if {$name ni $done} {
				putslog "name $name not found in $target yet ($testfile)"
				set ok 0
			}
		}
		if {!$ok} {
			putslog "$target has to be remade, renaming to $target.old"
			file rename -force -- $target $target.old
		}
	}
}
