proc cg_sam_clipamplicons {args} {
	if {([llength $args] != 3)} {
		errorformat sam_clipamplicons
		exit 1
	}
	foreach {ampliconsfile file outfile} {{} {} {}} break
	foreach {ampliconsfile file outfile} $args break
	set f [open $ampliconsfile]
	set header [tsv_open $f]
	close $f
	set poss [tsv_basicfields $header 3]
	lappend poss [lsearch $header outer_begin] [lsearch $header outer_end]
	if {[file extension $file] eq ".bam"} {
		set pipe [list samtools view -h $file]
	} else {
		set pipe [list {*}[gzcat $file] $file]
	}
	lappend pipe \| sam_clipamplicons $ampliconsfile {*}$poss
	if {[file extension $outfile] eq ".bam"} {
		lappend pipe \| samtools view -Shb - > $outfile
	} else {
		lappend pipe > $outfile
	}
	exec {*}$pipe
}
