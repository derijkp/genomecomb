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
