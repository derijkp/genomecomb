proc cg_regcommon {args} {
	if {([llength $args] != 1) && ([llength $args] != 2)} {
		errorformat regcommon
	}
	foreach {region_file1 region_file2} $args break
	set tempfile [tempfile]
	cg regsubtract $region_file1 $region_file2 > $tempfile
	cg regsubtract $region_file1 $tempfile >@ stdout
}
