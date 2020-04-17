proc samtools_sort {args} {
	if {[llength $args] < 2} {error  "format is: samtools_sort ?options? bamfile resultfile"}
	set bamfile [lindex $args end-1]
	set resultfile [lindex $args end]
	set args [lrange $args 0 end-2]
	if {[catch {version samtools 1}]} {
		if {[catch {exec samtools sort --no-PG {*}$args $bamfile $resultfile.temp 2>@ stdout} msg]} {
			error $msg
		}
		if {[file exists $resultfile.temp.bam]} {
			file rename -force -- $resultfile.temp.bam $resultfile
		} else {
			file rename -force -- $resultfile.temp $resultfile
		}
	} else {
		if {[catch {exec samtools sort --no-PG {*}$args $bamfile > $resultfile.temp 2>@ stdout} msg]} {
			error $msg
		}
		file rename -force -- $resultfile.temp $resultfile
	}
}
