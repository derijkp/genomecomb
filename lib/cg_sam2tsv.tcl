proc cg_sam2tsv {args} {
	set samfile -
	set outfile -
	set namecol name
	set informat sam
	cg_options sam2tsv args {
		-informat {set informat $value}
	} {samfile outfile} 0
	if {[file extension $samfile] in ".bam .cram"} {
		set informat bam
	} elseif {[file extension $samfile] eq ".sam"} {
		set informat sam
	}
	if {$informat eq "bam"} {
		if {$samfile eq "-"} {
			set pipe {samtools view -h | sam2tsv <@ stdin}
			set pipe 
		} else {
			set pipe [list samtools view -h $samfile \| sam2tsv]
		}
	} elseif {$samfile eq "-"} {
		set pipe {sam2tsv <@ stdin}
	} else {
		set pipe [list {*}[gzcat $samfile] $samfile \| sam2tsv]
	}
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

proc cg_bam2tsv {args} {
	cg_sam2tsv -informat bam {*}$args
}
