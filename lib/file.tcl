proc file_absolute {file} {
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

namespace eval genomecomb {}
if {[info commands genomecomb::cd.ori] eq ""} {
	rename cd genomecomb::cd.ori
}
if {[info commands genomecomb::pwd.ori] eq ""} {
	rename pwd genomecomb::pwd.ori
}
proc cd {path} {
	set ::genomecomb::cwd [file_absolute $path]
	genomecomb::cd.ori $::genomecomb::cwd
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
