proc cg_sam2tsv {args} {
	set samfile -
	set outfile -
	set namecol name
	cg_options sam2tsv args {
	} {samfile outfile} 0
	if {$samfile eq "-"} {
		set pipe {}
	} elseif {[file extension $samfile] eq ".bam"} {
		set pipe [list samtools view -h $samfile \|]
	} else {
		set pipe [list {*}[gzcat $samfile] $samfile \|]
	}
	lappend pipe sam2tsv
	if {$outfile eq "-"} {
		lappend pipe >@ stdout
	} else {
		lappend pipe > $outfile
	}
	if {$samfile eq "-"} {
		lappend pipe <@ stdin
	}
	catch_exec {*}$pipe
}
