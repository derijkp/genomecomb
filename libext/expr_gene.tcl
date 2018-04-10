proc tcl::mathfunc::codingcat {args} {
	if {[regexp CDS $args]} {
		return C
	} elseif {[regexp UTR $args]} {
		return U
	} elseif {[regexp RNA $args]} {
		return R
	} elseif {[regexp splice $args]} {
		return s
	} else {
		return -
	}
}

proc tcl::mathfunc::coding {args} {
	tcl::mathfunc::codingcat $args
}

proc tcl::mathfunc::zyg args {
	::zyg {*}$args
}

proc tcl::mathfunc::transcripts {genes impacts descrs filter format} {
	set result {}
	set descrs [split $descrs {,;}]
	set len [::llength $descrs]
	if {$len > 1} {
		set genes [split $genes {,;}]
		set impacts [split $impacts {,;}]
		if {[::llength $genes] == 1} {
			set genes [list_fill $len [lindex $genes 0]]
		}
		if {[::llength $impacts] == 1} {
			set impacts [list_fill $len [lindex $impacts 0]]
		}
		if {[::llength $descrs] > [::llength $genes]} {
			set genes [list_fill $len [lindex $genes 0]]
		}
	}
	foreach gene $genes descr $descrs impact $impacts {
		if {[::llength $filter] && $impact ni $filter} continue
		if {![regexp {^[+-]([^:]+):} $descr temp transcript]} {
			error "[lindex $args 0] has wrong format (should be a x_descr field)"
		}
		if {$gene eq "" || $format eq "t"} {
			lappend result $transcript
		} elseif {$format eq "g"} {
			lappend result $gene
		} elseif {$format eq "gt"} {
			lappend result $gene:$transcript
		} else {
			error "unknown format $format"
		}
	}
	return [join $result \;]
}

proc tcl::mathfunc::maximpact {args} {
	set list {}
	foreach el $args {
		lappend list {*}[split $el {,;}]
	}
	if {![::llength $list]} {return {}}
	set varlist [var_impact_list]
	set pos [::max [list_cor $varlist $list]]
	lindex $varlist $pos
}
