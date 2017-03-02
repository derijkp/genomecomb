proc samtools_sort {args} {
	if {[llength $args] < 2} {error  "format is: samtools_sort ?options? bamfile resultfile"}
	set bamfile [lindex $args end-1]
	set resultfile [lindex $args end]
	set args [lrange $args 0 end-2]
	if {[catch {version samtools 1}]} {
		exec samtools sort {*}$args $bamfile $resultfile 2>@ stderr
		file rename -force $resultfile.bam $resultfile
	} else {
		exec samtools sort {*}$args $bamfile > $resultfile 2>@ stderr
	}
}