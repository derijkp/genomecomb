proc cg_sam_clipamplicons {args} {
	if {([llength $args] != 3)} {
		errorformat sam_clipamplicons
	}
	foreach {ampliconsfile file outfile} {{} {} {}} break
	foreach {ampliconsfile file outfile} $args break
	set f [open $ampliconsfile]
	set header [tsv_open $f]
	close $f
	set poss [tsv_basicfields $header 3]
	lappend poss [lsearch $header outer_begin] [lsearch $header outer_end]
	if {[inlist $poss -1]} {
		error "error in amplicons file: missing fields: [list_sub {chromosome begin end outer_begin outer_end} [list_find $poss -1]]"
	}
	if {[file extension $file] eq ".bam"} {
		set pipe [list samtools view -h $file]
	} else {
		set pipe [list {*}[gzcat $file] $file]
	}
	lappend pipe \| sam_clipamplicons $ampliconsfile {*}$poss
	if {[file extension $outfile] eq ".bam" 
		|| ([file extension $outfile] eq ".temp" && [file extension [file root $outfile]] eq ".bam")} {
		lappend pipe \| samtools view -Shb - > $outfile
	} else {
		lappend pipe > $outfile
	}
	# putsvars pipe
	if {[catch {
		exec {*}$pipe
	} result]} {
		if {![regexp {^\[samopen\] SAM header is present: [0-9]+ sequences.$} $result]} {
			error $result
		}
	}
}
