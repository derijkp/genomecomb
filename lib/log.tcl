#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc logfile {logfile args} {
	set f [open $logfile a]
	foreach message $args {
		puts $f "[timestamp] $message"
	}
	close $f
}

set ::verbose 0

proc putslog {args} {
}

proc putsprogress {args} {
}

proc logverbose num {
	set ::verbose $num
	switch $num {
		0 {
			proc putslog {args} {
			}
			proc putsprogress {args} {
			}
			set ::cgjob(silent) 1
		}
		1 {
			proc putslog {args} {
				if {[lindex $args 0] eq "-nonewline"} {
					foreach message [lrange $args 1 end] {
						puts -nonewline stderr $message
					}
				} else {
					foreach message $args {
						puts stderr "[timestamp] $message"
					}
				}
			}
			proc putsprogress {args} {
			}
			set ::cgjob(silent) 0
		}
		2 {
			proc putslog {args} {
				if {[lindex $args 0] eq "-nonewline"} {
					foreach message [lrange $args 1 end] {
						puts -nonewline stderr $message
					}
				} else {
					foreach message $args {
						puts stderr "[timestamp] $message"
					}
				}
			}
			proc putsprogress {args} {
				if {[lindex $args 0] eq "-nonewline"} {
					foreach message [lrange $args 1 end] {
						puts -nonewline stderr $message
					}
				} else {
					foreach message $args {
						puts stderr "[timestamp] $message"
					}
				}
			}
			set ::cgjob(silent) 0
		}
	}
}
