proc tcl::mathfunc::simplifyschema {schema args} {
	set removeelements {}
	set trimelements {}
	if {[llength $args]} {
		foreach {removeelements trimelements} $args break
	}
	if {[llength $removeelements]} {set useremoveelements 1} else {set useremoveelements 0}
	if {[llength $trimelements]} {set usetrimelements 1} else {set usetrimelements 0}
	set simpleschema {}
	unset -nocomplain a
	set prev ""
	foreach {strand type} $schema {
		if {$useremoveelements && $type in $removeelements} continue
		if {$prev ne "" && $type eq $prev} continue
		incr a($strand)
		lappend simpleschema $type
		set prev $type
	}
	if {$usetrimelements} {
		set pos 0
		foreach el $simpleschema {
			if {$el ni $trimelements} break
			incr pos
		}
		set simpleschema [lrange $simpleschema $pos end]
		set pos [llength $simpleschema]
		foreach el [list_reverse $simpleschema] {
			incr pos -1
			if {$el ni $trimelements} break
		}
		set simpleschema [lrange $simpleschema 0 $pos]
	}
	if {[get a(-) 0] < [get a(+) 0]} {
		set simpleschema [list_reverse $simpleschema]
	}
	return $simpleschema
}
