proc info_analysis_file {resultfile sample parameters versions args} {
	set o [open $resultfile.temp w]
	puts $o [join {sample source parameter value} \t]
	# info
	foreach {item value} $args {
		puts $o [join [list $sample genomecomb $item $value] \t]
	}
	foreach param $parameters {
		puts $o [join [list $sample genomecomb param_$param [uplevel get $param ?]] \t]
	}
	foreach item $versions {
		puts $o [join [list $sample genomecomb version_$item [version $item]] \t]
	}
	close $o
	if {[file exists $resultfile]} {
		if {[catch {exec diff $resultfile $resultfile.temp}]} {
			file rename $resultfile $resultfile.[file_timestamp $resultfile]
		}
	}
	file rename -force -- $resultfile.temp $resultfile
}
