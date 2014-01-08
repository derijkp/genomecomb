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

proc tsv_select_spercent {arguments header neededfieldsVar} {
#	putsvars arguments
	upvar $neededfieldsVar neededfields
	foreach {cond1 cond2} $arguments break
	set samples [samples $header]
	set selection {}
	set total {}
	foreach sample $samples {
		set tempcond1 [tsv_select_replacevars $cond1 $header $sample]
		set tempcond2 [tsv_select_replacevars $cond2 $header $sample]
		lappend result \([tsv_select_detokenize $tempcond1 $header neededfields]\)
		lappend result \([tsv_select_detokenize $tempcond2 $header neededfields]\)
	}
	return "spercent_\([join $result ,]\)"
}

proc tsv_select_saggr {func1 func2 arguments header neededfieldsVar} {
	upvar $neededfieldsVar neededfields
# putsvars func1 func2 arguments header
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
	foreach sample $samples {
		if {$ifs ne ""} {
			set tempif [tsv_select_replacevars $ifs $header $sample]
			lappend result [tsv_select_detokenize $tempif $header neededfields]
		}
		set tempvalue [tsv_select_replacevars $values $header $sample]
		lappend result [tsv_select_detokenize $tempvalue $header neededfields]
	}
	return "$func\([join $result ,]\)"
}
