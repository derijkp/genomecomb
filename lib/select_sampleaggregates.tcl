proc tsv_select_replacevars {tokens header sample} {
	global newoptransa
	set result {}
	foreach line $tokens {
		foreach {type val} $line break
		switch -exact $type {
			@var {
				if {$val eq "sample"} {
					lappend result [list @val \"$sample\"]
				} elseif {[lsearch $header $val-$sample] != -1} {
					lappend result [list @var $val-$sample]
				} elseif {[lsearch $header $val] != -1} {
					lappend result [list @var $val]
				} else {
					lappend result [list @var $val-$sample]
				}
			}
			@function {
				set newline [lrange $line 0 1]
				set arguments [lrange $line 2 end]
				set ids {}
				foreach el $arguments {
					lappend newline [tsv_select_replacevars $el $header $sample]
				}
				lappend result $newline
			}
			@braces {
				set newline [lindex $line 0]
				lappend newline {*}[tsv_select_replacevars [lrange $line 1 end] $header $sample]
				lappend result $newline
			}
			default {
				lappend result $line
			}
		}
	}
	return $result
}

proc tsv_select_saggr_detokenize {tokens header neededfieldsVar missingVar} {
	upvar $neededfieldsVar neededfields
	upvar $missingVar missing
	if {[catch {tsv_select_detokenize $tokens $header neededfields} code]} {
		if {[regexp {field (.*) not present in file \(or sampleinfo\)} $code temp field]} {
			lappend missing $field
			return -code continue
		}
		error $code1
	}
	return $code
		
}

proc tsv_select_scount {arguments header neededfieldsVar} {
#	putsvars arguments
	upvar $neededfieldsVar neededfields
	set argument [lindex $arguments 0]
	set samples [samples $header]
	set missing {}
	set result {}
	foreach sample $samples {
		set temp [tsv_select_replacevars $argument $header $sample]
		set code [tsv_select_saggr_detokenize $temp $header neededfields missing]
		lappend result \($code\)
	}
	if {![llength $result]} {
		error "error in scount: all samples are missing one or more needed fields \([join $missing ,]\)"
	}
	return [join $result " + "]
}

proc tsv_select_spercent {arguments header neededfieldsVar} {
#	putsvars arguments
	upvar $neededfieldsVar neededfields
	foreach {cond1 cond2} $arguments break
	if {$cond1 eq ""} {error "error in spercent: at least one condition must be given (cannot be empty)"}
	if {$cond2 eq ""} {set cond2 $cond1 ; set cond1 1}
	set samples [samples $header]
	set selection {}
	set total {}
	set missing {}
	set result {}
	foreach sample $samples {
		if {$cond1 eq "1"} {
			set code1 1
		} else {
			set tempcond1 [tsv_select_replacevars $cond1 $header $sample]
			set code1 [tsv_select_saggr_detokenize $tempcond1 $header neededfields missing]
		}
		set tempcond2 [tsv_select_replacevars $cond2 $header $sample]
		set code2 [tsv_select_saggr_detokenize $tempcond2 $header neededfields missing]
		lappend result \($code1\)
		lappend result \($code2\)
	}
	if {![llength $result]} {
		error "error in spercent: all samples are missing one or more needed fields \([join $missing ,]\)"
	}
	return "spercent_\([join $result ,]\)"
}

proc tsv_select_saggr {func1 func2 arguments header neededfieldsVar} {
	# func1: function to be used if there is only 1 argument
	# func2: function to be used for 2 arguments (condition and value)
	upvar $neededfieldsVar neededfields
	if {[llength $arguments] == 1} {
		set values [lindex $arguments 0]
		set ifs {}
		set func $func1
	} else {
		foreach {ifs values} $arguments break
		set func $func2
	}
	set samples [samples $header]
	set result {}
	set missing {}
	foreach sample $samples {
		if {$ifs ne ""} {
			set tempif [tsv_select_replacevars $ifs $header $sample]
			set code1 [tsv_select_saggr_detokenize $tempif $header neededfields missing]
		}
		set tempvalue [tsv_select_replacevars $values $header $sample]
		set code2 [tsv_select_saggr_detokenize $tempvalue $header neededfields missing]
		if {$ifs ne ""} {
			lappend result $code1
		}
		lappend result $code2
	}
	if {![llength $result]} {
		error "error in [lindex [split $func _] 0]: all samples are missing one or more needed fields \([join $missing ,]\)"
	}
	return "$func\([join $result ,]\)"
}
