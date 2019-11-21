proc cg_sam_ampliconscount {args} {
	foreach {ampliconsfile file outfile} {{} {} {}} break
	cg_options sam_ampliconscount args {
	} {ampliconsfile file outfile} 2 3 {
		count reads mapping to each amplicon in ampliconsfile. file can be a bam or sam file.
	}
	set f [open $ampliconsfile]
	set header [tsv_open $f]
	close $f
	set poss [tsv_basicfields $header 3]
	lappend poss [lsearch $header outer_begin] [lsearch $header outer_end]
	if {[inlist $poss -1]} {
		error "error in amplicons file: missing fields: [list_sub {chromosome begin end outer_begin outer_end} [list_find $poss -1]]"
	}
	lappend poss [lsearch $header name]
	if {[file extension $file] in ".bam .cram"} {
		set pipe [list samtools view -h $file]
	} else {
		set pipe [list {*}[gzcat $file] $file]
	}
	lappend pipe \| sam_amplicons_count $ampliconsfile {*}$poss
	if {$outfile eq "" || $outfile eq "-"} {
		lappend pipe >@ stdout
	} else {
		lappend pipe {*}[compresspipe $outfile] > $outfile
	}
	# putsvars pipe
	if {[catch {
		 puts $pipe
		exec {*}$pipe
	} result]} {
		if {![regexp {^\[samopen\] SAM header is present: [0-9]+ sequences.$} $result]} {
			error $result
		}
	}
}
