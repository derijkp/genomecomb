proc tsv_select_replacevars {tokens header sample} {
	global newoptransa
	set result {}
	foreach line $tokens {
		foreach {type val} $line break
		switch -exact $type {
			@var {
				if {$val eq "sample"} {
					lappend result [list @val $sample]
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

proc tsv_select_scount {arguments header neededfieldsVar} {
#	putsvars arguments
	upvar $neededfieldsVar neededfields
	set argument [lindex $arguments 0]
	set samples [samples $header]
	set result {}
	foreach sample $samples {
		set temp [tsv_select_replacevars $argument $header $sample]
		lappend result \([tsv_select_detokenize $temp $header neededfields]\)
	}
	return [join $result " + "]
}

proc tsv_select_slist {arguments header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
	if {[llength $arguments] == 1} {
		set argument [lindex $arguments 0]
		set samples [samples $header]
		set result {}
		foreach sample $samples {
			set temp [tsv_select_replacevars $argument $header $sample]
			lappend result \[expr\ \{[tsv_select_detokenize $temp $header neededfields]\}\]
		}
		return \"[join $result ","]\"
	} else {
		foreach {ifs values} $arguments break
		set samples [samples $header]
		set result {}
		foreach sample $samples {
			set tempif [tsv_select_replacevars $ifs $header $sample]
			set tempvalue [tsv_select_replacevars $values $header $sample]
			lappend result [tsv_select_detokenize $tempif $header neededfields] \[expr\ \{[tsv_select_detokenize $tempvalue $header neededfields]\}\]
		}
		return "slist_if_\([join $result ,]\)"
	}
}
