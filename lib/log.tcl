#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc logfile {logfile args} {
	foreach message $args {
		set f [open $logfile a]
		puts $f "[timestamp] $message"
		close $f
	}
}

proc putslog {args} {
	foreach message $args {
		puts stderr "[timestamp] $message"
	}
}

proc putsprogress {args} {
	foreach message $args {
		puts stderr "[timestamp] $message"
	}
}

proc logverbose num {
	switch $num {
		0 {
			proc putslog {args} {
			}
			proc putsprogress {args} {
			}
		}
		1 {
			proc putslog {args} {
				foreach message $args {
					puts stderr "[timestamp] $message"
				}
			}
			proc putsprogress {args} {
			}
		}
		2 {
			proc putslog {args} {
				foreach message $args {
					puts stderr "[timestamp] $message"
				}
			}
			proc putsprogress {args} {
				foreach message $args {
					puts stderr "[timestamp] $message"
				}
			}
		}
	}
}
